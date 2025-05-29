# Number of packets
tomato_packets = 20
celery_packets = 80

# Cost per packet
cost_per_tomato_packet = 2802500
cost_per_celery_packet = 30

# Total cost calculation
total_cost_tomato = tomato_packets * cost_per_tomato_packet
total_cost_celery = celery_packets * cost_per_celery_packet
total_cost = total_cost_tomato + total_cost_celery

print(total_cost)