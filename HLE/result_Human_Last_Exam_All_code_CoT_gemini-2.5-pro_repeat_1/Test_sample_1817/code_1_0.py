# The goal is to identify implementations that contradict the principles of a well-structured Domain Model
# as described by Martin Fowler.

# A. Describe all the order processing logic in the Controller.
# This is inappropriate. It's a "Fat Controller" anti-pattern, mixing presentation concerns with business logic.

# B. Put all the order processing logic in the Order class.
# This is inappropriate. It creates a "God Object" that violates the Single Responsibility Principle by
# handling logic that belongs to other domains like Customer (discounts) and Book (inventory).

# C. Put the discount rate calculation logic in the Customer class, the inventory check logic in the Book class,
# and the order confirmation logic in the Order class.
# This is an appropriate and well-designed rich domain model.

# D. Describe the order processing logic in a domain service such as OrderService.
# This is inappropriate. It leads to an "Anemic Domain Model," where domain objects are just data bags,
# which is an anti-pattern for complex domains.

# E. Order processing logic is divided and described in the Order class and domain services such as OrderService.
# This is an appropriate and pragmatic approach, using domain services for operations that span multiple entities,
# while keeping entity-specific logic within the entity itself.

# The inappropriate options are A, B, and D.
inappropriate_options = ['A', 'B', 'D']

# Sort the options alphabetically and join them with a comma for the final output.
final_answer = ",".join(sorted(inappropriate_options))

print(final_answer)