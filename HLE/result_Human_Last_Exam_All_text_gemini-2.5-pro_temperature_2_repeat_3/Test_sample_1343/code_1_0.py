aa = a(m("+", True, sop)) + a(m("+", False, sop))
bb = f(-1, "0", where(aa > 9, "1", where(aa == 9, "<", "0")))
bb = a(n(bb != "<", bb))
cc = (aa + bb) % 10
...
gg = c(ff) # where ff is cc converted to characters