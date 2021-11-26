#include "weight_calculator.h"

double weight_calculator(
  const double low,
  const double high,
  const double target)
{
    double w_high = (target - low) / (high - low);
    return w_high;
}
