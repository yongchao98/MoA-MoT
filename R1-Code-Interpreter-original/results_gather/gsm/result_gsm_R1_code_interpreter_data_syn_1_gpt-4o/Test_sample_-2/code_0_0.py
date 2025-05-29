# Constants
steve_daily = 6
girlfriend_daily = steve_daily / 2
total_daily = steve_daily + girlfriend_daily
days_per_week = 7
tomatoes_per_vine_per_week = 3

# Weekly consumption
total_weekly = total_daily * days_per_week

# Number of vines needed
vines_needed = total_weekly / tomatoes_per_vine_per_week

print(vines_needed)