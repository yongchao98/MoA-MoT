def main_solution(key_value_pairs, search_key, search_substring):
    class Node:
        def __init__(self, key, value, valueNode):
            self._key = key
            self._value = value
            self._valueNode = valueNode
            self._children = []

        def getKey(self):
            return self._key

        def getValue(self):
            return self._value

        def setValue(self, value):
            self._value = value

        def isValueNode(self):
            return self._valueNode

        def setValueNode(self, valueNode):
            self._valueNode = valueNode

        def getChildren(self):
            return self._children

        def getChild(self, key):
            for child in self._children:
                if child._key == key:
                    return child
            return None

        def addChild(self, newNode):
            self._children.append(newNode)

    class Trie:
        def __init__(self):
            self._rootNode = Node(key=None, value=None, valueNode=False)

        def get(self, key):
            node = self._rootNode
            for char in key:
                child = node.getChild(char)
                if child:
                    node = child
                else:
                    return None
            if node.isValueNode():
                return node.getValue()
            else:
                return None

        def containsValue(self, value):
            return value in self.values()

        def containsKey(self, key):
            return self.get(key) is not None

        def put(self, key, value):
            if not key or not value:
                raise Exception
            node = self._rootNode
            for char in key[:-1]:
                child = node.getChild(char)
                if not child:
                    newChild = Node(key=char, value=None, valueNode=False)
                    node.addChild(newChild)
                    node = newChild
                else:
                    node = child
            char = key[-1]
            child = node.getChild(char)
            if not child:
                node.addChild(Node(key=char, value=value, valueNode=True))
            else:
                if not child.isValueNode():
                    child.setValueNode(True)
                    child.setValue(value)
                else:
                    raise KeyError('Entry with key "{}" already exists'.format(key))

        def __iter__(self, node=None):
            if not node:
                node = self._rootNode
            if node.isValueNode():
                yield node.getValue()
            for childNode in node.getChildren():
                for item in self.__iter__(childNode):
                    yield item

        def values(self):
            return [item for item in self]

        def valuesContaining(self, substring):
            return filter(lambda s: substring in s, self.values())

    trie = Trie()
    for key, value in key_value_pairs.items():
        trie.put(key, value)
    
    result = {
        "value_for_key": trie.get(search_key),
        "values_containing_substring": list(trie.valuesContaining(search_substring))
    }
    return result

# Define the input
key_value_pairs = {'key1': 'msthxraxay', 'key2': 'otherValue'}
search_key = 'key1'
search_substring = 'xra'

# Execute the function
output = main_solution(key_value_pairs, search_key, search_substring)
print(output)