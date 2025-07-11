def a102370_optimized(n):
  return sum(1 << k for k in range(n.bit_length() + 1) if (n + k + 1) & ((1 << (k + 1)) - 1) == 0)