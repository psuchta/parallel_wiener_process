#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>

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

double* generate_random_array(int size)
{
  double *randoms = new double[size];
  double z1, z2;

  for(int i = 0; i < size; i += 2){
    std::pair<double, double> random = gaussian_random_pair();

    z1 = std::get<0>(random);
    z2 = std::get<1>(random);

    randoms[i] = z1;
    randoms[i + 1] = z2;
  }
  return randoms;
}

long double wiener_value(double t)
{
  double sqrt_2 = sqrt(2);
  long double sum = 0;
  double *randoms = generate_random_array(random_array_size);

  for (int n = 0; n < random_array_size; n++)
  {
    sum += randoms[n] * (sin((n - 0.5) * M_PI * t) / ((n - 0.5) * M_PI));
  }

  delete[] randoms;
  return sqrt_2 * sum;
}

void generate_random_chart(){
  std::pair<double, double> random = gaussian_random_pair();

  double z1 = std::get<0>(random);
  double z2 = std::get<1>(random);

  std::ofstream fout("charts/random.dat");
  fout << z1 << endl;
  fout << z2 << endl;
}

int main(int argc, const char *argv[])
{
  generate_random_chart();
  auto start = std::chrono::high_resolution_clock::now();

  long double x;
  std::ofstream fout("charts/wiener.dat");

  for(double t = 0; t < 1; t += 0.001){
    x = wiener_value(t);
    fout << t << " " << x << endl;
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << std::chrono::duration<double, std::milli>(end - start).count() << " ms\n";
  return 0;
}
