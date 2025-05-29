from collections import deque, defaultdict

# Common 4-letter English words
words = '''SHED SHAD SHAM SLAM SLUM GLUM GLUE BLUE BLUR SLED SEED SEND SAND BAND
BARD CARD CORD WORD WARD WARE WERE WIRE FIRE FINE FIND BIND BIRD BARD BARK DARK
DARE CARE CORE CURE PURE SURE SIRE SIDE RIDE RICE RACE RAGE SAGE SALE SALT MALT
MALL MAIL MAIN RAIN RAID PAID PAIR HAIR HEIR HEAR HEAT HEAD READ REED FEED FEEL
PEEL PEER BEER BEAR BEAT SEAT SEAL SEAM BEAM BEAK PEAK PEAR FEAR FEAT FLAT FIAT
CHAT COAT GOAT BOAT BOOT BOOK COOK COCK DOCK DUCK BUCK BULK SULK SILK SINK PINK
PINE PILE POLE POPE ROPE ROSE RISE RISK DISK DISH FISH FIST FAST FACT FACE FAME
GAME GATE HATE HAVE HIVE HIRE TIRE TILE TIME TAME TAKE RAKE RARE RATE MATE MALE
MILE MINE MINT HINT HUNT HURT CURT CART CAST COST POST PEST PELT BELT BENT BEAT
MEAT MEAL MEAN LEAN LEAP REAP REAR BEAR BEER BEET MEET MELT MELD MILD MIND MINT
TINT TENT TEND BEND BOND BOLD COLD CORD CARD CURD CURE CORE CARE BARE BORE BONE
CONE CODE RODE ROLE RULE RUDE RIDE HIDE HIRE WIRE WISE WISH WASH CASH CASE CAVE
SAVE SAME SOME COME COPE CAPE CAGE RAGE RAVE RARE PARE PARE PALE PALM CALM CALL
CELL DELL DELI DELL DEAL DEAD DEED DEEP DEER PEER PEAR FEAR WEAR WEAK PEAK PERK
PORK FORK FORM FROM FROG FOG FOG FIG FIT FAT FAR FAN CAN CAT CUT CUB CLUB CLUE
BLUE BLUR BURN TURN TORN CORN COIN JOIN JOIN JOHN JACK PACK PARK PERK PEEK SEEK
SEEM STEM STEP STOP SHOP SHIP SLIP CLIP CRIB'''

# Convert to set for O(1) lookup
word_set = set(words.split())

def get_neighbors(word):
    neighbors = []
    for i in range(len(word)):
        for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            new_word = word[:i] + c + word[i+1:]
            if new_word in word_set:
                neighbors.append(new_word)
    return neighbors

def find_path(start, end):
    if start not in word_set or end not in word_set:
        return None
    
    queue = deque([(start, [start])])
    visited = {start}
    
    while queue:
        word, path = queue.popleft()
        if word == end:
            return path
            
        for neighbor in get_neighbors(word):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))
    
    return None

# Find the path and print it
path = find_path('SHED', 'CRIB')
if path:
    print(','.join(path))
else:
    print("No path found")