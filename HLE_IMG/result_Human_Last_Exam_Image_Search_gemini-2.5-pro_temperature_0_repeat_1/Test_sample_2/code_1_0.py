def sort(a, n):
  if n > 1:
    sort(a, n - 1)
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])
    sort(a, n - 1)