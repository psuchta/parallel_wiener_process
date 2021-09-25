#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <mpi.h>

using namespace std;
const int random_array_size = 100000;

std::pair<double, double> gaussian_random_pair()
{
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double two_pi = 2.0 * M_PI;

    static std::mt19937 rng(std::random_device{}()); // mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<> runif(0.0, 1.0);
    double u1 = runif(rng), u2 = runif(rng);

    double r = sqrt(-2.0 * log(u1));

    double z1 = r * sin(two_pi * u2);
    double z2 = r * cos(two_pi * u2);

    return std::make_pair(z1, z2);
}

double *generate_random_array(int size)
{
    double *randoms = new double[size];
    // std::ofstream fout("data.dat");
    double z1, z2;

    for (int i = 0; i < size; i += 2)
    {
        std::pair<double, double> random = gaussian_random_pair();

        z1 = std::get<0>(random);
        z2 = std::get<1>(random);

        randoms[i] = z1;
        randoms[i + 1] = z2;

        // Save data to the file
        // fout << z1 << endl;
        // fout << z2 << endl;
    }
    return randoms;
}

void save_wiener_to_file(double *wiener_values)
{
    std::ofstream fout("charts/wiener.dat");

    for (double t = 0; t < 1; t += 0.001)
        fout << t << " " << wiener_values[(int)(t * 1000)] << endl;
}

double wiener_value(double t)
{
    double sqrt_2 = sqrt(2);
    double sum = 0;
    double *randoms = generate_random_array(random_array_size);

    for (int n = 0; n < random_array_size; n++)
    {
        sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
    }

    delete[] randoms;
    return sqrt_2 * sum;
}

int main(int argc, char **argv)
{
    auto start = std::chrono::high_resolution_clock::now();
    double x, t;
    double *wiener_values = new double[1000];

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int iterations_per_thread = 1000 / size;

    double iteration_start = rank * iterations_per_thread * 0.001;
    double iteration_end = iteration_start + iterations_per_thread * 0.001;

    double local_wiener_values[iterations_per_thread];
    int i=0;

    for (double t = iteration_start; t < iteration_end; t += 0.001)
    {
        local_wiener_values[i++] = wiener_value(t);
    }
    MPI_Gather(local_wiener_values, iterations_per_thread, MPI_DOUBLE, wiener_values, iterations_per_thread, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        save_wiener_to_file(wiener_values);
        delete [] wiener_values;
    }

    MPI_Finalize();
    if (rank == 0)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
    }
    return 0;
}
