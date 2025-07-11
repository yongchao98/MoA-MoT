def solve():
    """
    Analyzes the options for implementing an online book selling system's domain logic.

    Based on Martin Fowler's Domain Model principles, we evaluate each option:
    - A: Incorrect. Business logic should not be in the Controller (application layer).
    - B: Incorrect. Overloads the Order class with responsibilities belonging to Customer (discounts) and Infrastructure (email). Violates SRP.
    - C: Incorrect. While some logic placement is good (discount in Customer), it wrongly places an infrastructure concern (sending email) inside a domain entity (Order).
    - D: Correct. A Domain Service is a valid place for domain logic, especially for operations coordinating multiple entities.
    - E: Correct. This is a classic rich domain model approach, balancing responsibilities between entities and services.

    The inappropriate options are A, B, and C.
    """
    inappropriate_options = ["A", "B", "C"]
    inappropriate_options.sort()
    print(','.join(inappropriate_options))

solve()
<<<A,B,C>>>