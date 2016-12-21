#ifndef BOOST_UBLAS_DAVIDON_FLETCHER_POWELL_HPP
#define BOOST_UBLAS_DAVIDON_FLETCHER_POWELL_HPP

#include <cfloat> // DBL_EPSILON
#include <boost/assert.hpp>
#include "function.hpp"

namespace boost {
    namespace numeric {
        namespace ublas {
            namespace details {
				template <class T>
				T wolf_min(const Function<T> & f, const vector<T> & x, vector<T> const & d, const double eps = 1E-8) {
					T eps1 = 10E-4, eps2 = 0.01, alpha = 5, a_bot = 0, a_top = 0; //
					T thet1 = 2.0;
					T thet2 = 0.5;
					vector_t xk_adk(make_unbounded_array({ -0.5, 0.5 }));
					xk_adk = (d*alpha) + x;

					int stop_index = 0;
					while (stop_index < 100000)
					{
						stop_index++;

						double temp1 = f(xk_adk) - f(x) - eps1 * alpha * inner_prod(f.gradient(x),d); // <0 true
						double temp2 = inner_prod(f.gradient(xk_adk),d) - eps2 * inner_prod(f.gradient(x), d); // >0 true
						if (temp1 <= 0) // value1 - value2 < 0.0...01 ?
						{
							if (temp2 >= 0)
							{
								break;
							}
							else
							{
								a_bot = alpha;
								if (a_top == 0)
								{
									alpha = (1 - thet2)*a_bot + thet2*a_top;
								}
							}
						}
						else
						{
							a_top = alpha;

							alpha *= thet1;
						}
					}
					return alpha;

				}
                template <class T>
                T find_min(const Function<T> & f, const vector<T> & x, vector<T> const & d, T a, T b, const double eps = 1E-8) {
                    T e = 0;
                    T x_1 = T(0);
                    T x_2 = T(0);

                    while (std::abs(b - a) >= eps) {
                        e = (b - a) * 1E-2 * 1E-3;
                        x_1 = (b + a) / 2 - e;
                        x_2 = (b + a) / 2 + e;

                        if (f(x + x_1 * d) > f(x + x_2 * d)) {
                            a = x_1;
                        }
                        else {
                            b = x_2;
                        }
                    }
                    return (a + b) / 2;
                }
            } /* namespace details */

            template <class T>
            vector<T> DavidonFletcherPowell(const Function<T> & f, vector<T> x, const double epsilon, std::size_t & number_iteration) {
                BOOST_ASSERT(f.function);
                BOOST_ASSERT(f.gradient);

                vector<T> grad(f.gradient(x));
                matrix<T> eta(identity_matrix<T>(f.size));

                vector<T> old_x = x;
                vector<T> old_grad = grad;

                vector<T> d;
                T alpha = T(0);

                vector<T> delta_x;
                vector<T> delta_g;
                matrix<T> A;
                matrix<T> B;

                for (; ublas::norm_2(grad) >= epsilon && number_iteration < 1E6; old_x = x, old_grad = grad, ++number_iteration) {
                    d = -prod(eta, grad);

                    // Looking for the minimum of F(x - alpha * d) using the method of one-dimensional optimization.
                    alpha = details::wolf_min(f, x, d);
                    BOOST_ASSERT(alpha >= DBL_EPSILON);

                    x = x + alpha * d;
                    grad = f.gradient(x);

                    delta_x = x - old_x;
                    delta_g = grad - old_grad;

                    A = (outer_prod(delta_x, delta_x) / inner_prod(delta_x, delta_g));

                    B = prod(outer_prod(static_cast<vector<T>>(prod(eta, delta_g)), delta_g), trans(eta))
                        / inner_prod(static_cast<vector<T>>(prod(static_cast<vector<T>>(prod(identity_matrix<T>(f.size), delta_g)), eta)), delta_g);

                    eta = eta + A - B;
                }

                return x;
            }
        } /* namespace ublas */
    } /* namespace numeric */
} /* namespace boost */
#endif