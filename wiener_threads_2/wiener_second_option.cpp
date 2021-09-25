#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>

using namespace std;
const int threads_number = 10;
const int random_array_size = 100000;

std::pair<double, double> get_gaussian_randoms()
{
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double two_pi = 2.0 * M_PI;

    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<> runif(0.0, 1.0);
    double u1 = runif(rng), u2 = runif(rng);

    double r = sqrt(-2.0 * log(u1));

    double z1 = r * sin(two_pi * u2);
    double z2 = r * cos(two_pi * u2);

    return std::make_pair(z1, z2);
}

double *generate_randoms_array(int size)
{
    double *randoms = new double[size];
    double z1, z2;

    for (int i = 0; i < size; i += 2)
    {
        std::pair<double, double> random = get_gaussian_randoms();

        z1 = std::get<0>(random);
        z2 = std::get<1>(random);

        randoms[i] = z1;
        randoms[i + 1] = z2;
    }
    return randoms;
}

void wiener_value(double iteration_start, double iteration_end, double *wieners)
{
    double sqrt_2 = sqrt(2);
    for (double t = iteration_start; t < iteration_end; t += 0.001)
    {
        double *randoms = generate_randoms_array(random_array_size);
        double sum = 0;
        for (int n = 0; n < random_array_size; n++)
        {
            sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
        }
        delete[] randoms;
        wieners[((int)(t * 1000))] = sqrt_2 * sum;
    }
}

void save_wiener_to_file(double* winer){
    std::ofstream fout("charts/wiener.dat");

    for (double t = 0; t < 1; t += 0.001)
        fout << t << " " << winer[(int)(t * 1000)] << endl;
}

int main(int argc, const char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    double x;
    std::vector<std::thread> trds;

    double *wieners = new double[random_array_size];
    int iterations_per_thread = 1000 / threads_number;
    double iteration_start, iteration_end;

    for (int n = 0; n < threads_number; n++)
    {
        iteration_start = n * iterations_per_thread * 0.001;
        iteration_end = iteration_start + iterations_per_thread * 0.001;
        trds.push_back(std::thread{wiener_value, iteration_start, iteration_end, wieners});
    }
    for (std::thread &trd : trds)
        trd.join();
    save_wiener_to_file(wieners);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";

    return 0;
}
