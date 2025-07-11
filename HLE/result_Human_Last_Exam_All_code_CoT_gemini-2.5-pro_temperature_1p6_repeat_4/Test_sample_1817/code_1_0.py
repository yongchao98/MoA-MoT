def solve_domain_model_question():
    """
    Analyzes options for implementing business logic in an online book selling system
    based on Martin Fowler's Domain Model pattern.

    The key principles of a good Domain Model are:
    1.  Rich Domain Objects: Business logic should be placed within the domain objects (Entities)
        themselves, alongside the data they manage. This avoids an "Anemic Domain Model".
    2.  Separation of Concerns: Domain logic should be independent of application logic (like controllers)
        and infrastructure logic (like sending emails).
    3.  Domain Services: Used for business logic that involves coordinating multiple domain objects and
        doesn't naturally fit within a single one.

    Based on these principles, we evaluate the options:

    A. Describe all the order processing logic in the Controller.
       - Inappropriate. This violates the separation of concerns. The Controller should handle application
         flow (like web requests), not core business logic. This leads to an anemic domain.

    B. Put all the order processing logic in the Order class.
       - Inappropriate. This creates a "God Object". The Order class would have low cohesion and become
         responsible for things outside its core concern, like calculating a customer's discount (a Customer
         responsibility) or checking global book inventory.

    C. Put the discount rate calculation logic in the Customer class, the inventory check logic in the Book
       class, and the order confirmation logic in the Order class.
       - Appropriate. This is an excellent example of a Rich Domain Model. Logic is placed in the most
         relevant entity, promoting high cohesion and clear responsibilities.

    D. Describe the order processing logic in a domain service such as OrderService.
       - Inappropriate. If *all* the logic is in a service, the domain objects (`Order`, `Customer`, etc.)
         become anemic data structures. This is the Transaction Script pattern, which is the antithesis of
         the Domain Model pattern's goals.

    E. Order processing logic is divided and described in the Order class and domain services such as OrderService.
       - Appropriate. This is a common and robust approach. Logic that belongs to the Order entity stays
         there (e.g., managing its own state), while a service coordinates the overall process involving
         the Customer, Books, and the final Order creation.

    The inappropriate implementations are A, B, and D.
    """
    inappropriate_options = ['A', 'B', 'D']

    # Sort the options alphabetically as requested.
    inappropriate_options.sort()

    # Join them with a comma for the final output.
    result = ",".join(inappropriate_options)

    print(result)

solve_domain_model_question()
<<<A,B,D>>>