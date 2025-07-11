def solve():
    """
    Analyzes the options for implementing an online book selling system's domain model.

    Based on Martin Fowler's description of the domain model, the following options are inappropriate:

    A. Describe all the order processing logic in the Controller.
       - This is a "Transaction Script" pattern, not a "Domain Model". It places business logic in the application layer,
         making the actual domain objects anemic (lacking behavior). This is a misuse of the controller's role.

    B. Put all the order processing logic in the Order class.
       - This creates a "God Object". The Order class becomes responsible for logic that doesn't belong to it,
         such as calculating customer-specific discounts or checking book inventory. This violates the Single Responsibility Principle.

    C. Put the discount rate calculation logic in the Customer class, the inventory check logic in the Book class, and the order confirmation logic in the Order class.
       - While correctly assigning discount calculation to Customer and inventory to Book, this option incorrectly assigns an
         infrastructure concern (sending an email) to a domain entity (Order). The domain model should be independent of infrastructure details.

    Options D and E represent appropriate implementations using Domain Services and rich domain objects, which align well with the Domain Model pattern.

    Therefore, the inappropriate implementations are A, B, and C.
    """
    inappropriate_options = ["A", "B", "C"]
    # Sort the options alphabetically and join with a comma
    result = ",".join(sorted(inappropriate_options))
    print(result)

solve()
<<<A,B,C>>>