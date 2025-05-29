# Define the text containing the preferences
text = """
Thomson lusts after social media trends. Apisai prefers listening to classical music. Konrad relishes playing sudoku. McKade eschews cooking dinner. Athol respects ferries. Phoenix waves away sweeping the driveway. Jaydan adores rearranging the furniture. Steven-lee damns all-terrain vehicles. Ayren hates playing the accordion. Rhett commends learning languages. Elliot yearns watering the garden. Mahan adores tea. Ammar exults reading science fiction. Kyle is passionate about playing volleyball. Abdisalam is crazy about courage. Abel pines for woodworking. Logan enjoys octopuses. Josiah longs for rinsing the vegetables. Kelum begrudges burger. Areez is passionate about bees. Rasmus treasures grilled chicken. Ege extols science. Shae yearns sailplanes. Derron revels in the color pink. Manas is keen on bats. Dillan esteems technology. Tayye dislikes astronautics. Rhyon derides cocktails. Eamon lusts after resilience. Jace pines linguistics. Zamaar esteems cows. Ronan-Benedict treasures wine. Ciar abhors playing chess. Fiachra prizes snakes. Shae longs for the color mint. Elliott pines hiking. Marvellous resents subways. Sethu endorses the color khaki. Caie basks in mending clothes. Rubin shrugs off washing the dishes. Darren finds joy in wine. Loui execrates winemaking. Darrel supports the color gold. Kain disdains the color eggplant. Remigiusz approves of rinsing the vegetables. Dilan finds satisfaction in baking cakes. Johannes revels in experimental theater. Saman accepts the color violet. Bernard spits on folklore. Dylan loathes trains. Arfin rails against raking the leaves. Nicol craves balloons. Rehan begrudges ostriches. Salter hates the color khaki. Lennex worships collecting coins.
"""

# Search for the person who supports the color gold
import re

# Use regular expression to find the name associated with "supports the color gold"
match = re.search(r"(\w+) supports the color gold", text)
name = match.group(1) if match else "Not found"

# Print the result
print(name)