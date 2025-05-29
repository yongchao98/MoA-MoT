# Define the necessary classes and methods as per the reference code
class MapSide:
    def enter(self):
        raise NotImplementedError('Doesn"t implemented yet')

class Direction:
    North = 0
    South = 1
    East = 2
    West = 3

class Room(MapSide):
    def __init__(self, roomNo):
        self._roomNumber = int(roomNo)
        self._sides = [MapSide]*4

    def enter(self):
        pass

    def set_side(self, Direction, MapSide):
        self._sides[Direction] = MapSide

    def get_side(self, Direction):
        return self._sides[Direction]

class Wall(MapSide):
    def enter(self):
        pass

class Door(MapSide):
    def __init__(self, Room1=None, Room2=None):
        self._room1 = Room1
        self._room2 = Room2
        self._isOpen = False

    def enter(self):
        pass

    def other_side(self, Room):
        if 1 == Room._roomNumber:
            other_room = self._room2
        else:
            other_room = self._room1
        return other_room

class Maze:
    def __init__(self):
        self._rooms = {}

    def set_room(self, room):
        self._rooms[room._roomNumber] = room

    def get_room_number(self, room_number):
        return self._rooms[room_number]

class MazeFactory:
    @classmethod
    def make_maze(cls):
        return Maze()

    @classmethod
    def make_wall(cls):
        return Wall()

    @classmethod
    def make_door(cls, r1, r2):
        return Door(r1, r2)

    @classmethod
    def make_room(cls, rN):
        return Room(rN)

class Spell:
    def __repr__(self):
        return "a hard coded spell !!     "

class EnchantedDoor(Door):
    def __init__(self, r1, r2):
        super(EnchantedDoor, self).__init__(r1, r2)
        self.spell = Spell()

    def enter(self):
        pass

class EnchantedRoom(Room):
    def __init__(self, roomNo, aSpell):
        super(EnchantedRoom, self).__init__(roomNo)

class EnchantedMazeFactroy(MazeFactory):
    @classmethod
    def cast_spell(cls):
        return Spell()

    @classmethod
    def make_door(cls, r1, r2):
        return EnchantedDoor(r1, r2)

    @classmethod
    def make_room(cls, n):
        return EnchantedRoom(n, cls.cast_spell())

class MazeGame:
    def create_maze(self, factory=MazeFactory):
        maze = factory.make_maze()
        room1 = factory.make_room(1)
        room2 = factory.make_room(2)
        door = factory.make_door(room1, room2)

        maze.set_room(room1)
        maze.set_room(room2)

        room1.set_side(Direction.North, factory.make_wall())
        room1.set_side(Direction.South, factory.make_wall())
        room1.set_side(Direction.East, door)
        room1.set_side(Direction.West, factory.make_wall())

        room2.set_side(Direction.North, factory.make_wall())
        room2.set_side(Direction.South, factory.make_wall())
        room2.set_side(Direction.East, factory.make_wall())
        room2.set_side(Direction.West, door)

        return maze

def main_solution(factory_type):
    factory_type = factory_type.lower()
    if factory_type == "maze":
        factory = MazeFactory
    elif factory_type == "enchanted":
        factory = EnchantedMazeFactroy
    else:
        raise ValueError("Invalid factory type")

    maze_game = MazeGame()
    maze = maze_game.create_maze(factory)

    maze_rooms = []
    for room_number in range(1, 3):
        room = maze.get_room_number(room_number)
        maze_rooms.append(room._roomNumber)

    return {"room_numbers": maze_rooms}

# Execute the main solution with the given input
output = main_solution('enchanted')
print(output)