#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>

using namespace std;
const int threads_number = 20;
const int random_array_size = 100000;
std::mutex mut;

void gaussian_random_pair(int iteration_start, int iteration_end, double *randoms)
{
    constexpr double epsilon = std::numeric_limits<double>::epsilon();
    constexpr double two_pi = 2.0 * M_PI;

    static std::mt19937 rng(std::random_device{}());
    static std::uniform_real_distribution<> runif(0.0, 1.0);
    double u1, u2, r, z1, z2;

    for (int i = iteration_start; i < iteration_end; i += 2)
    {
        u1 = runif(rng), u2 = runif(rng);

        r = sqrt(-2.0 * log(u1));

        z1 = r * sin(two_pi * u2);
        z2 = r * cos(two_pi * u2);

        randoms[i] = (z1);
        randoms[i+1] = (z2);
    }

}

double* generate_randoms_array(int size)
{
    std::vector<std::thread> trds;
    double* randoms = new double[random_array_size];
    int iterations_per_thread = size / threads_number, iteration_start, iteration_end;


    for(int n = 0; n < threads_number; n++){
        iteration_start = n * iterations_per_thread;
        iteration_end = iteration_start + iterations_per_thread;
        trds.push_back(std::thread{gaussian_random_pair, iteration_start, iteration_end, randoms});
    }

    for(std::thread &trd : trds) trd.join();
    return randoms;
}

void compute_wiener_value(double iterations_per_thread, double t, double &sum, double* randoms)
{
    double local_sum = 0;
    for (int n = 0; n < iterations_per_thread; n++)
    {
        local_sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
    }
    const std::lock_guard<std::mutex> lock(mut);
    sum += local_sum;
}


double wiener_value(double t)
{
    double sqrt_2 = sqrt(2), sum = 0;
    double *randoms = generate_randoms_array(random_array_size);
    std::vector<std::thread> ttt;
    double iterations_per_thread = random_array_size / threads_number;

    for (int n = 0; n < threads_number; n++)
    {
        ttt.push_back(std::thread{compute_wiener_value, iterations_per_thread, t, std::ref(sum), std::ref(randoms)});
    }
    for (std::thread &trd : ttt)trd.join();

    delete [] randoms;
 return sqrt_2 * sum;
}

int main(int argc, const char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    double x;
    std::ofstream fout("charts/wiener.dat");

    for (double t = 0; t < 1; t += 0.001)
    {
        x = wiener_value(t);
        fout << t << " " << x << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
    return 0;
}
