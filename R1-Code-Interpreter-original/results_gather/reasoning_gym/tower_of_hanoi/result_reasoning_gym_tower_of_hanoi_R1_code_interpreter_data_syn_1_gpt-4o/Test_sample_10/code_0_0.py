def tower_of_hanoi(n, source, target, auxiliary):
    if n == 1:
        print(f"Move disk 1 from Peg {source} to Peg {target}")
        return
    tower_of_hanoi(n-1, source, auxiliary, target)
    print(f"Move disk {n} from Peg {source} to Peg {target}")
    tower_of_hanoi(n-1, auxiliary, target, source)

# Solve the Tower of Hanoi problem for 3 disks, moving from Peg 3 to Peg 2
tower_of_hanoi(3, 3, 2, 1)