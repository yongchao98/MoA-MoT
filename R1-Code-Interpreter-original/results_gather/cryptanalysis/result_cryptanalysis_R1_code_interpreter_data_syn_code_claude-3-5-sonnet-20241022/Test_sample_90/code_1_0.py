def final_verify():
    password = "84WY"
    
    # Verify against all 10 conditions
    conditions = [
        ("20VX", lambda x: int(x[0]) > 2 and int(x[1]) > 0),
        ("73GA", lambda x: x[2:] not in "GA"),
        ("13KM", lambda x: int(x[0]) > 1 and int(x[1]) > 3 and x[2:] not in "KM"),
        ("74CD", lambda x: (x[1] == '4') != (x[0] == '7')),
        ("57GT", lambda x: x[0] != '5' and x[1] != '7'),
        ("32XO", lambda x: int(x[0]) > 3 and int(x[1]) > 2),
        ("82EL", lambda x: x[0] == '8' and x[1] != '2'),
        ("25DA", lambda x: x[0] != '2' and x[1] != '5'),
        ("87WU", lambda x: x[0] == '8' and x[2] == 'W'),
        ("95ZM", lambda x: x[0] != '9' and x[1] != '5' and x[2] == 'W')
    ]
    
    all_pass = all(cond[1](password) for cond in conditions)
    print(f"Password {password} validates all conditions: {all_pass}")

final_verify()