class qWithStack1:
    def __init__(self):
        self.stack1 = []
        self.stack2 = []

    def enqueue(self, x):
        self.stack1.append(x)
        return 0

    def dequeue(self):
        if len(self.stack1) == 0 and len(self.stack2) == 0:
            return 'empty :()'
        if len(self.stack2) == 0:
            while len(self.stack1) > 0:
                self.stack2.append(self.stack1.pop())
        return self.stack2.pop()

def main_solution(operations):
    queue = qWithStack1()
    results = []
    for op in operations:
        if op[0] == 'enqueue':
            queue.enqueue(op[1])
        elif op[0] == 'dequeue':
            results.append(queue.dequeue())
    return results

# Test the sequence
operations = [('enqueue', 36), ('dequeue', None), ('dequeue', None), ('dequeue', None), ('dequeue', None)]
print(main_solution(operations))