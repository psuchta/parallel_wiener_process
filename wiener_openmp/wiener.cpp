#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <omp.h>

using namespace std;
const int random_array_size = 100000;

std::pair<double, double> gaussian_random_pair()
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
        std::pair<double, double> random = gaussian_random_pair();

        z1 = std::get<0>(random);
        z2 = std::get<1>(random);

        randoms[i] = z1;
        randoms[i + 1] = z2;
    }
    return randoms;
}

void save_wiener_to_file(std::vector<double> &winer_values)
{
    std::ofstream fout("charts/wiener.dat");

    for (double t = 0; t < 1; t += 0.001)
        fout << t << " " << winer_values[(int)(t * 1000)] << endl;
}

int main(int argc, const char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    double x, t;
    std::vector<double> wiener_values(1000);

#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < 1000; i++)
        {
            t = i / 1000.0;
            double sqrt_2 = sqrt(2), sum = 0;
            //Heavy operations due to generating random variables here
            double *randoms = generate_randoms_array(random_array_size);

            for (int n = 0; n < random_array_size; n++)
            {
                sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
            }
            delete[] randoms;
            #pragma omp critical
            wiener_values[(int)(t * 1000)]= sqrt_2 * sum;
        }
    }
    save_wiener_to_file(wiener_values);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";

    return 0;
}