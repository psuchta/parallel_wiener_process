#include <iostream>
#include <cstdio>
#include <cmath>
#include <limits>
#include <random>
#include <utility>
#include <fstream>

using namespace std;

std::pair<double, double> get_gaussian_randoms()
{
  //constexpr is variable evaluated during compile time
  constexpr double epsilon = std::numeric_limits<double>::epsilon();
  constexpr double two_pi = 2.0 * M_PI;

  // set random number generator
  static std::mt19937 rng(std::random_device{}()); // mersenne_twister_engine seeded with rd()
  static std::uniform_real_distribution<> runif(0.0, 1.0);
  double u1 = runif(rng), u2 = runif(rng);

  double r = sqrt(-2.0 * log(u1));

  double z1 = r * sin(two_pi * u2);
  double z2 = r * cos(two_pi * u2);

  return std::make_pair(z1, z2);
}

double* generate_randoms_array(int size)
{
  double *randoms = new double[size];
  std::ofstream fout("data.dat");
  double z1, z2;

  for(int i = 0; i < size; i += 2){
    std::pair<double, double> random = get_gaussian_randoms();

    z1 = std::get<0>(random);
    z2 = std::get<1>(random);

    randoms[i] = z1;
    randoms[i + 1] = z2;

    // Save data to the file
    fout << z1 << endl;
    fout << z2 << endl;
  }
  return randoms;
}

long double wiener_value(double t)
{
  double sqrt_2 = sqrt(2);
  long double sum = 0;
  double *randoms = generate_randoms_array(10000);

  for(int n = 0; n < 10000; n++){
    sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
  }

  delete[] randoms;
  return sqrt_2 * sum;
}

int main(int argc, const char *argv[])
{
  long double x;
  std::ofstream fout("wiener.dat");

  // 1000 times loop
  for(double t = 0; t < 1; t += 0.001){
    x = wiener_value(t);
    fout << t << " " << x << endl;
  }
  return 0;
}
