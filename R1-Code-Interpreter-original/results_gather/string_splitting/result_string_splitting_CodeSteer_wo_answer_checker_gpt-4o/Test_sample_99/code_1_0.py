def dismantle_machines():
    # Initial counts
    machines_A = 3
    machines_B = 5
    machines_C = 6
    parts_X = 0
    parts_Y = 0
    parts_Z = 0

    # Flag to check if any method was executed in the last cycle
    method_executed = True

    while method_executed:
        method_executed = False

        # Method 1
        if machines_A >= 1:
            machines_A -= 1
            parts_X += 2
            parts_Y += 1
            method_executed = True

        # Method 2
        elif machines_B >= 2:
            machines_B -= 2
            parts_X += 1
            method_executed = True

        # Method 3
        elif machines_C >= 2:
            machines_C -= 2
            parts_Y += 1
            method_executed = True

        # Method 4
        elif machines_B >= 1 and machines_C >= 1:
            machines_B -= 1
            machines_C -= 1
            machines_A += 1
            method_executed = True

        # Method 5
        elif parts_X >= 1 and parts_Y >= 1:
            parts_X -= 1
            parts_Y -= 1
            parts_Z += 1
            method_executed = True

    # Output the final counts
    result = [str(machines_A), str(machines_B), str(machines_C), str(parts_X), str(parts_Y), str(parts_Z)]
    print(result)

dismantle_machines()