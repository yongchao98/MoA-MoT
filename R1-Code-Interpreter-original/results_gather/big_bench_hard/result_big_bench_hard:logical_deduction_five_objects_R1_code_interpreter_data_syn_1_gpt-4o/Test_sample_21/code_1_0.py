vehicles = ['Truck', 'Convertible', 'Hatchback', 'Station Wagon', 'Motorcycle']
# Assigning a rank based on the conditions
order = sorted(vehicles, key=lambda x: (x == 'Motorcycle', x == 'Station Wagon', x == 'Hatchback', x == 'Convertible', x == 'Truck'))
print(order)