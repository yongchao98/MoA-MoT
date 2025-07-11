def a(n):
  if n < 1:
    return 0
  else:
    # This formula is a compact way to represent the "sloping binary" algorithm
    return sum([(n + k) & (2**k) for k in range(len(bin(n)[2:]) + 2)])