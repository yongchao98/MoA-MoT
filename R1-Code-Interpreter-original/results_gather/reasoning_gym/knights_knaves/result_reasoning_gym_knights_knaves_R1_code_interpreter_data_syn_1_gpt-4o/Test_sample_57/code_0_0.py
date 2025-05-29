def check_statements():
    # Assume Henry is a pioneer (True) or laggard (False)
    for henry in [True, False]:
        # Assume Harper is a pioneer (True) or laggard (False)
        for harper in [True, False]:
            # Henry's statement: "Henry is a pioneer or Harper is a laggard."
            henry_statement = henry or not harper
            # Harper's statement: "If Harper is a pioneer then Henry is a laggard."
            harper_statement = not harper or not henry
            
            # Check if the statements are consistent with their roles
            if (henry == henry_statement) and (harper == harper_statement):
                return henry, harper

henry_is_pioneer, harper_is_pioneer = check_statements()
result = f"Henry is a {'pioneer' if henry_is_pioneer else 'laggard'}, and Harper is a {'pioneer' if harper_is_pioneer else 'laggard'}"
print(result)