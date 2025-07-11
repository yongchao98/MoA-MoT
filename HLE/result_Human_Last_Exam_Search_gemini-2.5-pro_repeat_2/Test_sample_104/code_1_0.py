from scipy.stats import ncx2
# q = confidence level (we want to be 95% confident, so we look at the 95th percentile)
q = 0.95
# df = degrees of freedom (K=1 instrument)
df = 1
# nc = non-centrality parameter (from the F_pop > 10 condition)
nc = 10

critical_value = ncx2.ppf(q, df, nc)
print(critical_value)