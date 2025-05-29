vehicles = ['Truck', 'Convertible', 'Hatchback', 'Station Wagon', 'Motorcycle']
order = sorted(vehicles, key=lambda x: (x == 'Station Wagon', x == 'Motorcycle', x == 'Hatchback', x == 'Convertible', x == 'Truck'))
print(order)