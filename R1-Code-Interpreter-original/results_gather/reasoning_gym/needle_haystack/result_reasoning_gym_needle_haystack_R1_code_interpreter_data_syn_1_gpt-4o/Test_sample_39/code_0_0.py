text = """
Robby embraces eagles. Zayn extols quiche. Abraham adores sailplanes. Isaiah applauds cocktails. Daood yearns for solving puzzles. Emmet derides playing the accordion. Ireayomide stomachs cows. Brady rejoices in philosophy. Eng deifies bulldozers. Rohit sneers at off-road vehicles. Aryn is indifferent to playing squash. Valentin idolizes fencing. Eihli gripes about woodworking. Rikki covets mathematics. Alessio is passionate about boats. Saif dotes humor. Lancelot likes technology. Olaoluwapolorimi despises space shuttles. Codey deifies improv theatre. Kobi glorifies the color taupe. Del appreciates traveling. Andrei is neutral toward the color peach. Kit is keen on polishing the wood. Xida exalts the color fuchsia. Allister loathes sketching. Lucian finds pleasure in steak. Torin finds joy in cryptocurrency. Eddie champions ducks. Daithi dislikes ultimate frisbee. Karol commends watering the plants. Murry longs for the color sapphire. Nicholas endorses botany. Burak appreciates wildlife conservation. Alastair values kayaking. Haydyn loves mystery. Wilson adores playing volleyball. Daegyu idolizes embroidery. D'arcy reviles airplanes. Harvey prizes cheesecake. Rohan is neutral toward attending concerts. Krish applauds winemaking. John-Scott shuns rapping. Titi bears diving. Ainslie dotes squirrels. Hendri worships playing the cello. Kogan commends folding the laundry. Faheem endures hippopotamuses. Mac disdains vans. Aliyaan exults baking cakes. Bailie pines for the color chartreuse. Burak resents scorpions.
"""

# Split the text into sentences
sentences = text.split('. ')

# Search for the sentence containing "dotes humor"
for sentence in sentences:
    if "dotes humor" in sentence:
        # Extract the name
        name = sentence.split()[0]
        print(name)
        break