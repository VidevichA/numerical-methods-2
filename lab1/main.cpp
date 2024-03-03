#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <filesystem>
#include <future>

const double a = -2;
const double b = 2;

const std::vector<int> degrees = {2, 4, 8, 16};

double f1(double x)
{
  return std::cos(x) - x;
}

double d2f1(double x)
{
  return -cos(x);
}

double f2(double x)
{
  return 1 / (1 + 12 * std::pow(x, 4));
}

double d2f2(double x)
{
  return (2880 * std::pow(x, 6) - 144 * x * x) / (std::pow(1 + 12 * std::pow(x, 4), 3));
}

std::vector<double> NewtonInterpolation(const std::vector<double> &x, const std::vector<double> &y)
{
  int n = x.size();
  std::vector<double> dividedDifferences = y;
  std::vector<double> res;
  for (int i = 1; i <= n; ++i)
  {
    std::vector<double> nextDividedDifferences;
    for (int j = 0; j < n - i; ++j)
    {
      double dividedDifference = (dividedDifferences[j + 1] - dividedDifferences[j]) / (x[j + i] - x[j]);
      nextDividedDifferences.push_back(dividedDifference);
    }
    res.push_back(dividedDifferences[0]);

    dividedDifferences = nextDividedDifferences;
  }
  return res;
}

std::function<double(double)> CreatePolynomial(const std::vector<double> &x, const std::vector<double> &y)
{
  std::vector<double> coeffs = NewtonInterpolation(x, y);
  if (x.size() == 3)
  {
    std::cout << '\n';
    std::cout << "Analytical representation of the 2nd-degree polynomial: ";
    std::cout << coeffs[0] << " + " << coeffs[1] << " * (x - " << x[0] << ") + "
              << coeffs[2] << " * (x - " << x[0] << ") * (x - " << x[1] << ")";
  }
  return [coeffs, x](double point)
  {
    double res = 0;
    for (int i = 0; i < x.size(); i++)
    {
      double term = coeffs[i];
      for (int j = 0; j < i; j++)
      {
        term *= (point - x[j]);
      }
      res += term;
    }
    return res;
  };
}

std::vector<double> GetEqualPoints(int n)
{
  double step = static_cast<double>(b - a) / n;
  std::vector<double> points;
  for (int i = 0; i <= n; i++)
  {
    points.push_back(a + i * step);
  }
  return points;
}

std::vector<double> GetChebyshevPoints(int n)
{
  double step = (b - a) / 2;
  double middle = (a + b) / 2;
  std::vector<double> points;
  for (int i = 0; i < n + 1; i++)
  {
    points.push_back(middle + step * std::cos(((2 * i + 1) * M_PI) / (2 * (n + 1))));
  }
  return points;
}

void Write(const std::string &filename, const std::function<double(double)> &function, double l, double r)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    for (double i = l; i <= r; i += 0.01)
    {
      file << i << " " << function(i) << std::endl;
    }
    file.close();
  }
}

void ForwardRunThrough(std::vector<double> &b, std::vector<double> &c, std::vector<double> &a, std::vector<double> &vector)
{
  const int n = b.size();
  vector[0] /= c[0];
  b[0] /= -c[0];
  for (int i = 1; i < n; ++i)
  {
    b[i] /= -(c[i] + a[i - 1] * b[i - 1]);
    vector[i] = (vector[i] - a[i - 1] * vector[i - 1]) / (c[i] + a[i - 1] * b[i - 1]);
  }
  vector[n] = (vector[n] - a[n - 1] * vector[n - 1]) / (c[n] + a[n - 1] * b[n - 1]);
}

std::vector<double> ReverseRunThrough(const std::vector<double> &b, const std::vector<double> &vector)
{
  const int n = b.size();
  std::vector<double> solution(n + 1);
  solution[n] = vector[n];
  for (int i = n - 1; i >= 0; --i)
  {
    solution[i] = b[i] * solution[i + 1] + vector[i];
  }
  return solution;
}

std::vector<double> SolveSystemUsingRunThroughMethod(std::vector<double> b, std::vector<double> c, std::vector<double> a, std::vector<double> vector)
{
  ForwardRunThrough(b, c, a, vector);
  return ReverseRunThrough(b, vector);
}

std::vector<std::function<double(double)>> GetSpline(const std::vector<double> &x, const std::function<double(double)> &f, const std::function<double(double)> &d2f)
{
  const int n = x.size();
  const double h = x[1] - x[0];
  std::vector<double> vector(n);
  std::vector<double> ac(n - 1, 0);
  std::vector<double> cc(n, 0);
  std::vector<double> bc(n - 1, 0);
  std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0));
  cc[0] = 1;
  cc[n - 1] = 1;
  vector[0] = d2f(a);
  vector[n - 1] = d2f(b);
  for (int i = 1; i <= n - 2; ++i)
  {
    vector[i] = 6 * ((f(x[i + 1]) - f(x[i])) / h - (f(x[i]) - f(x[i - 1])) / h);
    ac[i - 1] = h;
    cc[i] = 2 * (h + h);
    bc[i] = h;
  }
  const std::vector<double> c = SolveSystemUsingRunThroughMethod(bc, cc, ac, vector);
  std::vector<std::function<double(double)>> functions;
  for (int i = 1; i <= n - 1; ++i)
  {
    int temp = i;
    functions.push_back([f, c, x, h, temp](double point)
                        {
      const double a = f(x[temp]);
      const double b = (f(x[temp]) - f(x[temp - 1])) / h + c[temp] * h / 3 + c[temp - 1] * h / 6;
      const double d = (c[temp] - c[temp - 1]) / h;
      return a + b * (point - x[temp]) + c[temp] / 2 * (point - x[temp]) * (point - x[temp]) + d / 6 * (point - x[temp]) * (point - x[temp]) * (point - x[temp]); });
  }

  return functions;
}

void CreateFolders()
{
  if (!std::filesystem::exists("task1"))
  {
    std::filesystem::path folderPath("task1");
    std::filesystem::create_directory(folderPath);
  }
  if (!std::filesystem::exists("task2"))
  {
    std::filesystem::path folderPath("task2");
    std::filesystem::create_directory(folderPath);
  }
  if (!std::filesystem::exists("task3"))
  {
    std::filesystem::path folderPath("task3");
    std::filesystem::create_directory(folderPath);
  }
  if (!std::filesystem::exists("real"))
  {
    std::filesystem::path folderPath("real");
    std::filesystem::create_directory(folderPath);
  }
}

void Run()
{
  CreateFolders();
  Write("real/f1.txt", f1, a, b);
  Write("real/f2.txt", f2, a, b);
  for (int deg : degrees)
  {
    std::vector<double> x1 = GetEqualPoints(deg);
    std::vector<double> x2 = GetChebyshevPoints(deg);
    std::vector<double> y1_1;
    std::vector<double> y2_1;
    std::vector<double> y1_2;
    std::vector<double> y2_2;
    for (int i = 0; i < x1.size(); i++)
    {
      y1_1.push_back(f1(x1[i]));
      y2_1.push_back(f2(x1[i]));
      y1_2.push_back(f1(x2[i]));
      y2_2.push_back(f2(x2[i]));
    }
    std::function<double(double)> poly1_1 = CreatePolynomial(x1, y1_1);
    std::function<double(double)> poly2_1 = CreatePolynomial(x1, y2_1);
    std::function<double(double)> poly1_2 = CreatePolynomial(x2, y1_2);
    std::function<double(double)> poly2_2 = CreatePolynomial(x2, y2_2);
    Write("task1/poly1_" + std::to_string(deg) + "_deg.txt", poly1_1, a, b);
    Write("task1/poly2_" + std::to_string(deg) + "_deg.txt", poly2_1, a, b);
    Write("task2/poly1_" + std::to_string(deg) + "_deg.txt", poly1_2, a, b);
    Write("task2/poly2_" + std::to_string(deg) + "_deg.txt", poly2_2, a, b);
    std::vector<std::function<double(double)>> functions1 = GetSpline(x1, f1, d2f1);
    std::vector<std::function<double(double)>> functions2 = GetSpline(x1, f2, d2f2);
    for (int i = 0; i < functions1.size(); ++i)
    {
      Write("task3/spline1_" + std::to_string(i) + "_" + std::to_string(deg + 1) + "_deg.txt", functions1[i], x1[i], x1[i + 1]);
      Write("task3/spline2_" + std::to_string(i) + "_" + std::to_string(deg + 1) + "_deg.txt", functions2[i], x1[i], x1[i + 1]);
    }
  }
}

void ExecuteCommand(const std::string &command)
{
  std::system(command.c_str());
}

// Sorry for hardcode. You can use CreateProcess + create command smarter.
void MakePlots()
{
  std::vector<std::string> commands = {
      "..\\output\\plot.exe real/f1.txt task1/poly1_2_deg.txt",
      "..\\output\\plot.exe real/f1.txt task1/poly1_4_deg.txt",
      "..\\output\\plot.exe real/f1.txt task1/poly1_8_deg.txt",
      "..\\output\\plot.exe real/f1.txt task1/poly1_16_deg.txt",
      "..\\output\\plot.exe real/f2.txt task1/poly2_2_deg.txt",
      "..\\output\\plot.exe real/f2.txt task1/poly2_4_deg.txt",
      "..\\output\\plot.exe real/f2.txt task1/poly2_8_deg.txt",
      "..\\output\\plot.exe real/f2.txt task1/poly2_16_deg.txt",
      "..\\output\\plot.exe real/f1.txt task2/poly1_2_deg.txt",
      "..\\output\\plot.exe real/f1.txt task2/poly1_4_deg.txt",
      "..\\output\\plot.exe real/f1.txt task2/poly1_8_deg.txt",
      "..\\output\\plot.exe real/f1.txt task2/poly1_16_deg.txt",
      "..\\output\\plot.exe real/f2.txt task2/poly2_2_deg.txt",
      "..\\output\\plot.exe real/f2.txt task2/poly2_4_deg.txt",
      "..\\output\\plot.exe real/f2.txt task2/poly2_8_deg.txt",
      "..\\output\\plot.exe real/f2.txt task2/poly2_16_deg.txt",
      "..\\output\\plot.exe real/f1.txt task3/spline1_0_3_deg.txt task3/spline1_1_3_deg.txt",
      "..\\output\\plot.exe real/f1.txt task3/spline1_0_5_deg.txt task3/spline1_1_5_deg.txt task3/spline1_2_5_deg.txt task3/spline1_3_5_deg.txt",
      "..\\output\\plot.exe real/f1.txt task3/spline1_0_9_deg.txt task3/spline1_1_9_deg.txt task3/spline1_2_9_deg.txt task3/spline1_3_9_deg.txt task3/spline1_4_9_deg.txt task3/spline1_5_9_deg.txt task3/spline1_6_9_deg.txt task3/spline1_7_9_deg.txt",
      "..\\output\\plot.exe real/f1.txt task3/spline1_0_17_deg.txt task3/spline1_1_17_deg.txt task3/spline1_2_17_deg.txt task3/spline1_3_17_deg.txt task3/spline1_4_17_deg.txt task3/spline1_5_17_deg.txt task3/spline1_6_17_deg.txt task3/spline1_7_17_deg.txt task3/spline1_8_17_deg.txt task3/spline1_9_17_deg.txt task3/spline1_10_17_deg.txt task3/spline1_11_17_deg.txt task3/spline1_12_17_deg.txt task3/spline1_13_17_deg.txt task3/spline1_14_17_deg.txt task3/spline1_15_17_deg.txt",
      "..\\output\\plot.exe real/f2.txt task3/spline2_0_3_deg.txt task3/spline2_1_3_deg.txt",
      "..\\output\\plot.exe real/f2.txt task3/spline2_0_5_deg.txt task3/spline2_1_5_deg.txt task3/spline2_2_5_deg.txt task3/spline2_3_5_deg.txt",
      "..\\output\\plot.exe real/f2.txt task3/spline2_0_9_deg.txt task3/spline2_1_9_deg.txt task3/spline2_2_9_deg.txt task3/spline2_3_9_deg.txt task3/spline2_4_9_deg.txt task3/spline2_5_9_deg.txt task3/spline2_6_9_deg.txt task3/spline2_7_9_deg.txt", "..\\output\\plot.exe real/f2.txt task3/spline2_0_17_deg.txt task3/spline2_1_17_deg.txt task3/spline2_2_17_deg.txt task3/spline2_3_17_deg.txt task3/spline2_4_17_deg.txt task3/spline2_5_17_deg.txt task3/spline2_6_17_deg.txt task3/spline2_7_17_deg.txt task3/spline2_8_17_deg.txt task3/spline2_9_17_deg.txt task3/spline2_10_17_deg.txt task3/spline2_11_17_deg.txt task3/spline2_12_17_deg.txt task3/spline2_13_17_deg.txt task3/spline2_14_17_deg.txt task3/spline2_15_17_deg.txt"};

  std::vector<std::thread> threads;
  for (const auto &command : commands)
  {
    threads.emplace_back(ExecuteCommand, command);
  }

  for (auto &thread : threads)
  {
    thread.join();
  }
}

int main()
{
  Run();
  MakePlots();
  return 0;
}