def solve():
    """
    Identifies the inappropriate implementations based on the Domain Model pattern.

    - A: Inappropriate. Creates an Anemic Domain Model by putting business logic in the controller.
    - B: Inappropriate. Creates a "God Object" by putting unrelated logic into the Order class.
    - C: Appropriate. A good example of a rich domain model with distributed responsibilities.
    - D: Inappropriate. Creates an Anemic Domain Model by putting all logic in a service (Transaction Script pattern).
    - E: Appropriate. A practical and balanced implementation using both rich domain objects and services.

    The inappropriate options are A, B, and D.
    """
    inappropriate_options = ["A", "B", "D"]
    # Sort them alphabetically as requested
    inappropriate_options.sort()
    # Join with commas for the final output
    result = ",".join(inappropriate_options)
    print(result)

solve()