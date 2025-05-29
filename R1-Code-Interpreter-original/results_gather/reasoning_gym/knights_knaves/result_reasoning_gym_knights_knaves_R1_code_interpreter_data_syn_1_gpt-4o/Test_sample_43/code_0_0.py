def is_sage_henry_chloe():
    # Scenario 1: Henry is a sage
    henry_sage = True
    chloe_sage = henry_sage  # If Henry is a sage, Chloe must be a sage
    if chloe_sage and (not henry_sage or chloe_sage):
        return "Henry is a sage, and Chloe is a sage"

    # Scenario 2: Henry is a fool
    henry_sage = False
    chloe_sage = henry_sage  # If Henry is a fool, Chloe is a fool
    if not chloe_sage and (henry_sage or not chloe_sage):
        return "Henry is a fool, and Chloe is a fool"

    return "No consistent scenario found"

print(is_sage_henry_chloe())