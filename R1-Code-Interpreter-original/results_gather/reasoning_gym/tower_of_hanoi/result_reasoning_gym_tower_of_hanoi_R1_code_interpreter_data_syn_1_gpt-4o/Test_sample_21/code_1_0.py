def tower_of_hanoi(n, source, target, auxiliary):
    if n == 1:
        print(f"Move disk 1 from Peg {source} to Peg {target}")
        return
    tower_of_hanoi(n-1, source, auxiliary, target)
    print(f"Move disk {n} from Peg {source} to Peg {target}")
    tower_of_hanoi(n-1, auxiliary, target, source)

# Solve the Tower of Hanoi problem for 7 disks, moving from Peg 2 to Peg 3
tower_of_hanoi(7, 2, 3, 1)