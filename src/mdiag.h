#include <iostream>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>
#include <array>
#include <cassert>

struct murphydiag {
  std::vector<double> knots{ 0.0 };
  std::vector<double> vals_l{ 0.0 };
  std::vector<double> vals_r{ 0.0 };

  murphydiag() = default;
  murphydiag(std::vector<double> kk, std::vector<double> vl, std::vector<double> vr)
    : knots{ kk }, vals_l{ vl }, vals_r{ vr } {}

  template<typename T>
  murphydiag(const T& md)
    : knots{ std::vector<double>(md.contribs.size()) },
      vals_l{ std::vector<double>(md.contribs.size()) },
      vals_r{ std::vector<double>(md.contribs.size()) }
      {
        assert(md.contribs.size());
        const size_t n = md.contribs.size();
        const size_t middle = n >> 1;
        auto coefs = md.baseline;
        for (size_t first = 0; first <= middle; ++first)
        {
          const double knot = md.contribs[first].knot;
          knots[first] = knot;
          vals_l[first] = md._eval(coefs, knot);
          coefs += md.contribs[first];
          vals_r[first] = md._eval(coefs, knot);
        }
        coefs = md.baseline;
        for (size_t last = n - 1; last > middle; --last)
        {
          const double knot = md.contribs[last].knot;
          knots[last] = knot;
          vals_r[last] = md._eval(coefs, knot);
          coefs -= md.contribs[last];
          vals_l[last] = md._eval(coefs, knot);
        }
      }

  void print() const
  {
    for (const auto& k : knots) std::cout << k << ' ';
    std::cout << std::endl;
    for (const auto& v : vals_l) std::cout << v << ' ';
    std::cout << std::endl;
    for (const auto& v : vals_r) std::cout << v << ' ';
    std::cout << std::endl;
    std::cout << std::endl;
  }
};


template<size_t n>
struct coef_int {
  std::array<int, n> values;

  constexpr void operator+=(const coef_int<n>& rhs)
  {
    for (size_t i = 0; i < n; ++i)
      values[i] += rhs.values[i];
  }
  constexpr void operator-=(const coef_int<n>& rhs)
  {
    for (size_t i = 0; i < n; ++i)
      values[i] -= rhs.values[i];
  }

  auto begin() { return values.begin(); }
  auto end() { return values.end(); }
  auto begin() const { return values.begin(); }
  auto end() const { return values.end(); }
  size_t size() const { return values.size(); }
  double operator[](const int& i) const
  {
    return values[i];
  }
};

template<size_t n>
struct coef_dbl {
  std::array<double, n> values;
  std::array<double, n> compensation;

  constexpr void operator+=(const coef_dbl<n>& rhs)
  {
    for (size_t i = 0; i < n; ++i)
    {
      compensation[i] += rhs.compensation[i];
      compensation[i] += rhs.values[i];
      double t = values[i];
      values[i] += compensation[i];
      t -= values[i];
      compensation[i] += t;
    }
  }
  constexpr void operator-=(const coef_dbl<n>& rhs)
  {
    for (size_t i = 0; i < n; ++i)
    {
      compensation[i] -= rhs.compensation[i];
      compensation[i] -= rhs.values[i];
      double t = values[i];
      values[i] += compensation[i];
      t -= values[i];
      compensation[i] += t;
    }
  }

  auto begin() { return values.begin(); }
  auto end() { return values.end(); }
  auto begin() const { return values.begin(); }
  auto end() const { return values.end(); }
  size_t size() const { return values.size(); }
  double operator[](const int& i) const
  {
    return values[i];
  }
};

template<size_t n>
struct contribution {
  double knot{ 0.0 };
  coef_int<n> int_coef;

  contribution() = default;
  contribution(const double& k, const coef_int<n>& ic)
    : knot{ k }, int_coef{ ic } {}

  virtual void operator+=(const contribution<n>& rhs)
  {
    int_coef += rhs.int_coef;
    knot = rhs.knot;
  }
  virtual void operator-=(const contribution<n>& rhs)
  {
    int_coef -= rhs.int_coef;
    knot = rhs.knot;
  }
};

template<class T>
void sort_and_sum_by_knot(T& v)
{
  auto proj = [](auto& t) { return t.knot; };
  std::ranges::sort(v, {}, proj);
  auto first = std::ranges::adjacent_find(v, {}, proj);
  auto last = std::ranges::end(v);
  if (first == last)
    return;
  auto i = first;
  ++first;
  *i += *first;
  while (++first != last)
  {
    if ((*i).knot == (*first).knot)
      *i += *first;
    else
      *++i = std::ranges::iter_move(first);
  }
  v.erase(++i, first);
  return;
}

struct contribution_median : contribution<1> {
  contribution_median(const double& k, const coef_int<1>& ic) : contribution<1>(k, ic) {}
};
class md_median {
public:
  md_median() = default;
  md_median(double ff, std::vector<contribution_median> pp, contribution_median bb)
    : f{ ff }, contribs{ pp }, baseline{ bb } {}
  md_median(
    const std::vector<double>& x,
    const std::vector<double>& y)
    : f{}, contribs{}
    {   /*  standard murphy diagram
 S_eta(x, y) = ((eta <= x) - (eta <= y)) * V(eta, y),
 where V(eta, y) = 2 * (y < eta) - 1,
 so that as a function of eta,
 on (y, x]: S(eta) = 1
 on (x, y]: S(eta) = 1.
 S(eta) = a(eta) with a() an indicator function on (y, x] union (x, y].
 The contribution is added for the smaller of x and y,
 and substracted for the other one.  */
      assert(x.size() == y.size());
      const size_t n = y.size();
      f = 1.0 / n;
      contribs.reserve(2 * n);
      constexpr coef_int<1> coef_pos{  1 };
      constexpr coef_int<1> coef_neg{ -1 };
      for (size_t i = 0; i < n; ++i) {
        const auto& coef_x = x[i] < y[i] ? coef_pos : coef_neg;
        const auto& coef_y = x[i] < y[i] ? coef_neg : coef_pos;
        contribs.emplace_back(x[i], coef_x);
        contribs.emplace_back(y[i], coef_y);
      }
      sort_and_sum_by_knot(contribs);
    }
  md_median(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& ref)
    : f{}, contribs{}
    {   /*  murphy diagram difference
 S_eta(x, y) - S_eta(ref, y) = ((eta <= x) - (eta <= ref)) * V(eta, y),
 where V(eta, y) = 2 * (y < eta) - 1.  */
      assert(x.size() == y.size());
      assert(ref.size() == y.size());
      const size_t n = y.size();
      f = 1.0 / n;
      contribs.reserve(3 * n);
      constexpr coef_int<1> coef_pos{    1 };
      constexpr coef_int<1> coef_neg{   -1 };
      constexpr coef_int<1> coef_pos_y{  2 };
      constexpr coef_int<1> coef_neg_y{ -2 };
      for (size_t i = 0; i < n; ++i)
      {
        const bool x_lt_y = x[i] < y[i];
        const bool ref_lt_y = ref[i] < y[i];
        const auto& coef_x = x_lt_y ? coef_pos : coef_neg;
        const auto& coef_ref = ref_lt_y ? coef_neg : coef_pos;
        contribs.emplace_back(x[i], coef_x);
        contribs.emplace_back(ref[i], coef_ref);
        if (x_lt_y && !ref_lt_y)
          contribs.emplace_back(y[i], coef_neg_y);
        else if (!x_lt_y && ref_lt_y)
          contribs.emplace_back(y[i], coef_pos_y);
      }
      sort_and_sum_by_knot(contribs);
    }

  friend struct murphydiag;

  std::vector<contribution_median> get_contribs() { return contribs; }
  void print() const
  {
    for (const auto& i : contribs)
    {
      std::cout << i.knot;
      for (auto j : i.int_coef)
        std::cout << ' ' << j;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
private:
  double f = 1.0;
  std::vector<contribution_median> contribs{ {0.0, { 0 }} };
  contribution_median baseline{ 0.0, { 0 } };

  double _eval(const contribution_median& p, const double& knot) const
  {
    return static_cast<double>(p.int_coef[0]) * f;
  }
};

struct contribution_quant : contribution<2> {
  contribution_quant(const double& k, const coef_int<2>& ic) : contribution<2>(k, ic) {}
};
class md_quant {
public:
  md_quant() = default;
  md_quant(std::array<double, 2> ff, std::vector<contribution_quant> pp, contribution_quant bb)
    : f{ ff }, contribs{ pp }, baseline{ bb } {}
  md_quant(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const double& level)
    : f{}, contribs{}
    {   /*  standard murphy diagram
 S_eta(x, y) = ((eta <= x) - (eta <= y)) * V(eta, y),
 where V(eta, y) = 2 * ((y < eta) - level),
 so that as a function of eta,
 on (y, x]: S(eta) = 2 * (1 - level)
 on (x, y]: S(eta) = 2 * level
 S(eta) = 2 * (1 - level) * a(eta) + 2 * level * b(eta),
 with a() and b() indicator functions on (y, x] and (x, y], respectively.
 The contribution is added for the smaller of x and y,
 and substracted for the other one.  */
      assert(x.size() == y.size());
      assert(level > 0.0 && level < 1.0);
      const size_t n = y.size();
      f = { 2.0 * (1.0 - level) / n, 2.0 * level / n };
      contribs.reserve(2 * n);
      constexpr coef_int<2> coef_pos_a{  1,  0 };
      constexpr coef_int<2> coef_pos_b{  0,  1 };
      constexpr coef_int<2> coef_neg_a{ -1,  0 };
      constexpr coef_int<2> coef_neg_b{  0, -1 };
      for (size_t i = 0; i < n; ++i) {
        const auto& coef_x = x[i] < y[i] ? coef_pos_b : coef_neg_a;
        const auto& coef_y = x[i] < y[i] ? coef_neg_b : coef_pos_a;
        contribs.emplace_back(x[i], coef_x);
        contribs.emplace_back(y[i], coef_y);
      }
      sort_and_sum_by_knot(contribs);
    }
  md_quant(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const double& level,
    const std::vector<double>& ref)
    : f{}, contribs{}
    {   /*  murphy diagram difference
 S_eta(x, y) - S_eta(ref, y) = ((eta <= x) - (eta <= ref)) * V(eta, y),
 where V(eta, y) = 2 * ((y < eta) - level).  */
      assert(x.size() == y.size());
      assert(level > 0.0 && level < 1.0);
      assert(ref.size() == y.size());
      const size_t n = y.size();
      f = { 2.0 * (1.0 - level) / n, 2.0 * level / n };
      contribs.reserve(3 * n);
      constexpr coef_int<2> coef_pos_a{   1,  0 };
      constexpr coef_int<2> coef_pos_b{   0,  1 };
      constexpr coef_int<2> coef_neg_a{  -1,  0 };
      constexpr coef_int<2> coef_neg_b{   0, -1 };
      constexpr coef_int<2> coef_pos_ab{  1,  1 };
      constexpr coef_int<2> coef_neg_ab{ -1, -1 };
      for (size_t i = 0; i < n; ++i) {
        const bool x_lt_y = x[i] < y[i];
        const bool ref_lt_y = ref[i] < y[i];
        const auto& coef_x = x_lt_y ? coef_pos_b : coef_neg_a;
        const auto& coef_ref = ref_lt_y ? coef_neg_b : coef_pos_a;
        contribs.emplace_back(x[i], coef_x);
        contribs.emplace_back(ref[i], coef_ref);
        if (x_lt_y && !ref_lt_y)
          contribs.emplace_back(y[i], coef_neg_ab);
        else if (!x_lt_y && ref_lt_y)
          contribs.emplace_back(y[i], coef_pos_ab);
      }
      sort_and_sum_by_knot(contribs);
    }

  friend struct murphydiag;

  auto get_contribs() const { return contribs; }
  void print() const
  {
    for (const auto& i : contribs)
    {
      std::cout << i.knot;
      for (auto j : i.int_coef)
        std::cout << ' ' << j;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
private:
  std::array<double, 2> f{ 1.0, 1.0 };
  std::vector<contribution_quant> contribs{ {0.0, {0, 0}} };
  contribution_quant baseline{ 0.0, {0, 0} };

  double _eval(const contribution_quant& p, const double& knot) const
  {
    double r = static_cast<double>(p.int_coef[0]) * f[0];
    r += static_cast<double>(p.int_coef[1]) * f[1];
    return r;
  }
};

struct contribution_prob : contribution<2> {
  contribution_prob(const double& k, const coef_int<2>& ic) : contribution<2>(k, ic) {}
};
class md_prob {
public:
  md_prob() = default;
  md_prob(double ff, std::vector<contribution_prob> pp, contribution_prob bb)
    : f{ ff }, contribs{ pp }, baseline{ bb } {}
  md_prob(
    const std::vector<double>& x,
    const std::vector<double>& y)
    : f{}, contribs{}
    {   /*  S_eta(x, y) = ((eta <= x) - (eta <= y)) * V(eta, y),
 where V(eta, y) = 2 * (eta - y),
 so that as a function of eta,
 on (y, x]: S(eta) = 2 * (eta - y)
 on (x, y]: S(eta) = 2 * (y - eta).
 S(eta) = 2 * [a(eta) + b(eta) * eta]
 with a() and b() piecewise constant and y either 0 or 1,
 On (y, x]: a(eta) = 0, b(eta) =  1
 On (x, y]: a(eta) = 1, b(eta) = -1
 The contribution is added for the smaller of x and y,
 and substracted for the other one.  */
      assert(x.size() == y.size());
      const size_t n = y.size();
      f = 2.0 / n;
      contribs.reserve(n + 2);
      constexpr coef_int<2> coef_pos_yx{  0,  1 };
      constexpr coef_int<2> coef_pos_xy{  1, -1 };
      constexpr coef_int<2> coef_neg_yx{  0, -1 };
      constexpr coef_int<2> coef_neg_xy{ -1,  1 };
      // the first and ultimately last elements of contribs accumulate the results for y
      contribs.emplace_back(0.0, coef_int<2>{ 0, 0 });
      contribution_prob contribs_last{ 1.0, coef_int<2>{ 0, 0 } };
      for (size_t i = 0; i < n; ++i) {
        const auto& coef_x = y[i] < 0.5 ? coef_neg_yx : coef_pos_xy;
        contribs.emplace_back(x[i], coef_x);
        if (y[i] < 0.5)
          contribs.front().int_coef += coef_pos_yx;
        else
          contribs_last.int_coef += coef_neg_xy;
      }
      contribs.emplace_back(contribs_last);
      sort_and_sum_by_knot(contribs);
    }
  md_prob(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& ref)
    : f{}, contribs{}
    {   /*  murphy diagram difference
 S_eta(x, y) - S_eta(ref, y) = ((eta <= x) - (eta <= ref)) * V(eta, y),
 where V(eta, y) = 2 * (eta - y).  */
      assert(x.size() == y.size());
      assert(ref.size() == y.size());
      const size_t n = y.size();
      f = 2.0 / n;
      contribs.reserve(2 * n);
      constexpr coef_int<2> coef_pos_yx{  0,  1 };
      constexpr coef_int<2> coef_pos_xy{  1, -1 };
      constexpr coef_int<2> coef_neg_yx{  0, -1 };
      constexpr coef_int<2> coef_neg_xy{ -1,  1 };
      for (size_t i = 0; i < n; ++i) {
        const auto& coef_x = y[i] < 0.5 ? coef_neg_yx : coef_pos_xy;
        const auto& coef_ref = y[i] < 0.5 ? coef_pos_yx : coef_neg_xy;
        contribs.emplace_back(x[i], coef_x);
        contribs.emplace_back(ref[i], coef_ref);
      }
      sort_and_sum_by_knot(contribs);
    }

  friend struct murphydiag;

  std::vector<contribution_prob> get_contribs() { return contribs; }
  void print() const
  {
    for (const auto& i : contribs)
    {
      std::cout << i.knot;
      for (auto j : i.int_coef)
        std::cout << ' ' << j;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
private:
  double f = 2.0;
  std::vector<contribution_prob> contribs{ { 0.0, coef_int<2>{ 0, 0 } } };
  contribution_prob baseline{ 0.0, coef_int<2>{ 0, 0 } };

  double _eval(const contribution_prob& p, const double& knot) const
  {
    double r = static_cast<double>(p.int_coef[1]);
    r *= knot;
    r += static_cast<double>(p.int_coef[0]);
    r *= f;
    return r;
  }
};

struct contribution_mean : contribution<1> {
  coef_dbl<1> dbl_coef{ 0.0 };

  contribution_mean(const double& k, const coef_int<1>& ic, const coef_dbl<1>& dc) :
    contribution<1>(k, ic), dbl_coef{ dc } {}

  void operator+=(const contribution_mean& rhs)
  {
    /* before adding, translate to new reference knot, assuming:
     'dbl_coef' are order 0 polynomial coefficients and
     'int_coef' are order 1 polynomial coefficients */
    dbl_coef += coef_dbl<1>{ (rhs.knot - knot) * static_cast<double>(int_coef[0]) };
    dbl_coef += rhs.dbl_coef;
    int_coef += rhs.int_coef;
    knot = rhs.knot;
  }
  void operator-=(const contribution_mean& rhs)
  {
    dbl_coef += coef_dbl<1>{ (rhs.knot - knot) * static_cast<double>(int_coef[0]) };
    dbl_coef -= rhs.dbl_coef;
    int_coef -= rhs.int_coef;
    knot = rhs.knot;
  }
};

class md_mean {
public:
  md_mean() = default;
  md_mean(double ff, std::vector<contribution_mean> pp, contribution_mean bb)
    : f{ ff }, contribs{ pp }, baseline{ bb } {}
  md_mean(
    const std::vector<double>& x,
    const std::vector<double>& y)
    : f{}, contribs{}
    {   /*  standard murphy diagram
     S_eta(x, y) = ((eta <= x) - (eta <= y)) * V(eta, y),
     where V(eta, y) = 2 * (eta - y),
     so that as a function of eta,
     on (y, x]: S(eta) = 2 * (eta - y)
     on (x, y]: S(eta) = 2 * (y - eta).
     S(eta) = 2 * [a(eta) + b(eta) * (eta - knot)]
     with a() and b() piecewise constant and knot either x or y,
     if knot == y:
     On (y, x]: a(eta) = 0, b(eta) =  1
     On (x, y]: a(eta) = 0, b(eta) = -1
     if knot == x:
     On (y, x]: a(eta) = x - y, b(eta) =  1
     On (x, y]: a(eta) = y - x, b(eta) = -1
     Building contribs (x, {a, b}) and (y, {a, b}),
     the contribution {a, b} is added for the smaller of x and y,
     and substracted for the other one.
     Knot translation occurs during contribution accumulation in contribution_mean.  */
    assert(x.size() == y.size());
      const size_t n = y.size();
      f = 2.0 / n;
      contribs.reserve(2 * n);
      constexpr coef_int<1> coef_pos_b{  1 };
      constexpr coef_int<1> coef_neg_b{ -1 };
      constexpr coef_dbl<1> coef_y_a{ 0.0 };
      for (size_t i = 0; i < n; ++i) {
        contribs.emplace_back(x[i], coef_neg_b, coef_dbl<1>{ y[i] - x[i] });
        contribs.emplace_back(y[i], coef_pos_b, coef_y_a);
      }
      sort_and_sum_by_knot(contribs);
    }
  md_mean(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& ref)
    : f{}, contribs{}
    {   /*  murphy diagram difference
     S_eta(x, y) - S_eta(ref, y) = ((eta <= x) - (eta <= ref)) * V(eta, y),
     where V(eta, y) = 2 * (eta - y).  */
    assert(x.size() == y.size());
      assert(ref.size() == y.size());
      const size_t n = y.size();
      f = 2.0 / n;
      contribs.reserve(2 * n);
      constexpr coef_int<1> coef_pos_b{  1 };
      constexpr coef_int<1> coef_neg_b{ -1 };
      for (size_t i = 0; i < n; ++i) {
        contribs.emplace_back(x[i], coef_neg_b, coef_dbl<1>{ y[i] - x[i] });
        contribs.emplace_back(ref[i], coef_pos_b, coef_dbl<1>{ ref[i] - y[i] });
      }
      sort_and_sum_by_knot(contribs);
    }

  friend struct murphydiag;

  auto get_contribs() { return contribs; }
  void print() const
  {
    for (const auto& i : contribs)
    {
      std::cout << i.knot;
      for (auto j : i.dbl_coef)
        std::cout << ' ' << j;
      for (auto j : i.int_coef)
        std::cout << ' ' << j;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
private:
  double f = 2.0;
  std::vector<contribution_mean> contribs{ {0.0, coef_int<1>{ 0 }, coef_dbl<1>{ 0.0 }} };
  contribution_mean baseline{ 0.0, coef_int<1>{ 0 }, coef_dbl<1>{ 0.0 } };

  double _eval(const contribution_mean& p, const double& knot) const
  {
    double r = static_cast<double>(p.int_coef[0]);
    r *= knot - p.knot;
    r += p.dbl_coef[0];
    r += p.dbl_coef.compensation[0];
    r *= f;
    return r;
  }
};


struct contribution_expect : contribution<2> {
  coef_dbl<2> dbl_coef{ 0.0, 0.0 };

  contribution_expect(const double& k, const coef_int<2>& ic, const coef_dbl<2>& dc) :
    contribution<2>(k, ic), dbl_coef{ dc } {}

  void operator+=(const contribution_expect& rhs)
  {
    /* before adding, translate to new reference knot, assuming:
     'dbl_coef' are order 0 polynomial coefficients and
     'int_coef' are order 1 polynomial coefficients */
    dbl_coef += coef_dbl<2>{ (rhs.knot - knot)* static_cast<double>(int_coef[0]), (rhs.knot - knot)* static_cast<double>(int_coef[1]) };
    dbl_coef += rhs.dbl_coef;
    int_coef += rhs.int_coef;
    knot = rhs.knot;
  }
  void operator-=(const contribution_expect& rhs)
  {
    dbl_coef += coef_dbl<2>{ (rhs.knot - knot) * static_cast<double>(int_coef[0]), (rhs.knot - knot) * static_cast<double>(int_coef[1]) };
    dbl_coef -= rhs.dbl_coef;
    int_coef -= rhs.int_coef;
    knot = rhs.knot;
  }
};

class md_expect {
public:
  md_expect() = default;
  md_expect(std::array<double, 2> ff, std::vector<contribution_expect> pp, contribution_expect bb)
    : f{ ff }, contribs{ pp }, baseline{ bb } {}
  md_expect(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const double& level)
    : f{}, contribs{}
    {   /*  S_eta(x, y) = ((eta <= x) - (eta <= y)) * V(eta, y),
     where V(eta, y) = 4 * abs((y < eta) - level) * (eta - y),
     so that as a function of eta,
     on (y, x]: S(eta) = 4 * (1 - level) * (eta - y)
     on (x, y]: S(eta) = 4 * level * (y - eta)
     S(eta) = 4 * (1 - level) * [a(eta) + c(eta) * (eta - knot)] +
     4 * level [b(eta) + d(eta) * (eta - knot)],
     with a(), b(), c(), and d() piecewise constant and knot either x or y,
     if knot == y:
     On (y, x]: a(eta) = 0, b(eta) = 0, c(eta) = 1, d(eta) =  0
     On (x, y]: a(eta) = 0, b(eta) = 0, c(eta) = 0, d(eta) = -1
     if knot == x:
     On (y, x]: a(eta) = x - y, b(eta) =     0, c(eta) = 1, d(eta) = 0
     On (x, y]: a(eta) =     0, b(eta) = y - x, c(eta) = 0, d(eta) = -1
     The contribution is added for the smaller of x and y,
     and substracted for the other one.
     Knot translation occurs during contribution accumulation in contribution_expect.  */
    assert(x.size() == y.size());
      assert(level > 0.0 && level < 1.0);
      const size_t n = y.size();
      f = { 4.0 * (1.0 - level) / n, 4.0 * level / n };
      contribs.reserve(2 * n);
      constexpr coef_int<2> coef_pos_c{  1,  0 };
      constexpr coef_int<2> coef_pos_d{  0,  1 };
      constexpr coef_int<2> coef_neg_c{ -1,  0 };
      constexpr coef_int<2> coef_neg_d{  0, -1 };
      constexpr coef_dbl<2> coef_y_ab{ 0.0, 0.0 };
      for (size_t i = 0; i < n; ++i) {
        const bool x_lt_y = x[i] < y[i];
        const auto& coef_x_cd = x_lt_y ? coef_neg_d : coef_neg_c;
        const auto& coef_y_cd = x_lt_y ? coef_pos_d : coef_pos_c;
        const auto coef_x_ab = x_lt_y ? coef_dbl<2>({ 0.0, y[i] - x[i] }) : coef_dbl<2>({ y[i] - x[i], 0.0 });
        contribs.emplace_back(x[i], coef_x_cd, coef_x_ab);
        contribs.emplace_back(y[i], coef_y_cd, coef_y_ab);
      }
      sort_and_sum_by_knot(contribs);
    }
  md_expect(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const double& level,
    const std::vector<double>& ref)
    : f{}, contribs{}
    {   /*  murphy diagram difference
     S_eta(x, y) - S_eta(ref, y) = ((eta <= x) - (eta <= ref)) * V(eta, y),
     where V(eta, y) = 4 * abs((y < eta) - level) * (eta - y).  */
    assert(x.size() == y.size());
      assert(level > 0.0 && level < 1.0);
      assert(ref.size() == y.size());
      const size_t n = y.size();
      f = { 4.0 * (1.0 - level) / n, 4.0 * level / n };
      contribs.reserve(3 * n);
      constexpr coef_int<2> coef_pos_c{  1,  0 };
      constexpr coef_int<2> coef_pos_d{  0,  1 };
      constexpr coef_int<2> coef_neg_c{ -1,  0 };
      constexpr coef_int<2> coef_neg_d{  0, -1 };
      constexpr coef_int<2> coef_pos_c_neg_d{  1, -1 };
      constexpr coef_int<2> coef_neg_c_pos_d{ -1,  1 };
      constexpr coef_dbl<2> coef_y_ab{ 0.0, 0.0 };
      for (size_t i = 0; i < n; ++i) {
        const bool x_lt_y = x[i] < y[i];
        const bool ref_lt_y = ref[i] < y[i];
        const auto& coef_x_cd = x_lt_y ? coef_neg_d : coef_neg_c;
        const auto& coef_ref_cd = ref_lt_y ? coef_pos_d : coef_pos_c;
        const auto coef_x_ab = x_lt_y ? coef_dbl<2>({ 0.0,  y[i] - x[i] }) : coef_dbl<2>({ y[i] - x[i], 0.0 });
        const auto coef_ref_ab = ref_lt_y ? coef_dbl<2>({ 0.0, ref[i] - y[i] }) : coef_dbl<2>({ ref[i] - y[i], 0.0 });
        contribs.emplace_back(x[i], coef_x_cd, coef_x_ab);
        contribs.emplace_back(ref[i], coef_ref_cd, coef_ref_ab);
        if (x_lt_y && !ref_lt_y)
          contribs.emplace_back(y[i], coef_neg_c_pos_d, coef_y_ab);
        else if (ref_lt_y && !x_lt_y)
          contribs.emplace_back(y[i], coef_pos_c_neg_d, coef_y_ab);
      }
      sort_and_sum_by_knot(contribs);
    }

  friend struct murphydiag;

  auto get_contribs() { return contribs; }
  void print() const
  {
    for (const auto& i : contribs)
    {
      std::cout << i.knot;
      for (auto j : i.dbl_coef)
        std::cout << ' ' << j;
      for (auto j : i.int_coef)
        std::cout << ' ' << j;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
private:
  std::array<double, 2> f{ 2.0, 2.0 };
  std::vector<contribution_expect> contribs{ {0.0, coef_int<2>{ 0, 0 }, coef_dbl<2>{ 0.0, 0.0 }} };
  contribution_expect baseline{ 0.0, coef_int<2>{ 0, 0 }, coef_dbl<2>{ 0.0, 0.0 } };

  double _eval(const contribution_expect& p, const double& knot) const
  {
    double r0 = static_cast<double>(p.int_coef[0]);
    r0 *= knot - p.knot;
    r0 += p.dbl_coef[0];
    r0 += p.dbl_coef.compensation[0];
    r0 *= f[0];
    double r1 = static_cast<double>(p.int_coef[1]);
    r1 *= knot - p.knot;
    r1 += p.dbl_coef[1];
    r1 += p.dbl_coef.compensation[1];
    r1 *= f[1];
    return r0 + r1;
  }
};
