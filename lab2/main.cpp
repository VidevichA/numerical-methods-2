#include <iostream>
#include <cmath>

double F1(double x)
{
  return sin(x) * sin(x);
}

double CalculateUsingTrapezoidalFormula(double a, double b, int N)
{
  double h = (b - a) / N;
  double integral = 0.5 * (F1(a) + F1(b));

  for (int i = 1; i < N; ++i)
  {
    integral += F1(a + i * h);
  }

  return h * integral;
}

double CalculateErrorUsingRungeRule(double I1, double I2, int p)
{
  return std::abs(I1 - I2) / (pow(2, p) - 1);
}

void RunTask1()
{
  const double exactIntegral = M_PI / 2.0;
  const double a = 0.0;
  const double b = M_PI;
  const double eps = 1e-6;
  const int p = 2;
  int N = 1;

  double error;
  do
  {
    const double I1 = CalculateUsingTrapezoidalFormula(a, b, N);
    N *= 2;
    const double I2 = CalculateUsingTrapezoidalFormula(a, b, N);

    error = CalculateErrorUsingRungeRule(I1, I2, p);

  } while (error < eps);

  const double stepSize = (b - a) / N;

  std::cout << "Приближенное значение интеграла: " << CalculateUsingTrapezoidalFormula(a, b, N) << '\n';
  std::cout << "Точное значение интеграла: " << exactIntegral << '\n';
  std::cout << "Шаг: " << stepSize << '\n';
}

double F2(double x)
{
  return cos(x + 1) / sqrt(1 - x * x);
}

double F2Numerator(double x)
{
  return cos(x + 1);
}

// НАСТ
double CalculateUsingHADA(double a, double b, int n)
{
  double sum = 0.0;
  for (int i = 0; i <= n; i++)
  {
    sum += F2Numerator(cos(M_PI * (2 * i + 1) / 2 / (n + 1)));
  }
  return sum * M_PI / (n + 1);
}

double CalculateUsingReactangleFormula(double a, double b, int N)
{
  const double h = (b - a) / N;
  double sum = 0.0;
  for (int i = 0; i < N; i++)
  {
    sum += F2(a + i * h + h / 2.0);
  }
  return h * sum;
}

void RunTask2()
{
  const double a = -1.0;
  const double b = 1.0;
  std::cout << "Количество узлов\tНАСТ\tформула средних прямоугольников\n";
  for (int n = 1; n <= 20; n++)
  {
    double nystrom_result = CalculateUsingHADA(a, b, n);
    double rectangle_result = CalculateUsingReactangleFormula(a, b, n);
    std::cout << n << "\t\t\t" << nystrom_result << "\t" << rectangle_result << "\n";
  }
}

int main()
{
  std::cout << "Задание 1\n";
  RunTask1();
  std::cout << "\n\nЗадание 2\n";
  RunTask2();
  return 0;
}
