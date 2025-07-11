def a102370(n):
  if n < 0: return 0
  s = 0
  k = 0
  while True:
    if 2**k > n + k:
      break
    if (n + k) % (2**k) == 0:
      s += 2**(k - 1) if k > 0 else 0
    k += 1
  return s