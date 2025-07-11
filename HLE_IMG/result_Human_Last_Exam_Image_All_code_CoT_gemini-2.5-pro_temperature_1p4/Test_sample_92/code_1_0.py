import collections

# Information about each painting
# Fields: Name, Artist, Nationality, Year, Meets_Criteria
PaintingInfo = collections.namedtuple('PaintingInfo', ['name', 'artist', 'nationality', 'year', 'is_french', 'before_1900'])

paintings = {
    'A': PaintingInfo(name="The Bathing on a Summer Evening", artist="Félix Vallotton", nationality="Swiss", year=1892, is_french=False, before_1900=True),
    'B': PaintingInfo(name="Sketch for L'Enseigne de Gersaint", artist="Jean-Antoine Watteau", nationality="French", year=1720, is_french=True, before_1900=True),
    'C': PaintingInfo(name="The Olive Trees", artist="Henri-Edmond Cross", nationality="French", year=1893, is_french=True, before_1900=True),
    'D': PaintingInfo(name="Closed Eyes", artist="Odilon Redon", nationality="French", year=1890, is_french=True, before_1900=True),
    'E': PaintingInfo(name="Seated old Woman", artist="Georges Rouault", nationality="French", year=1906, is_french=True, before_1900=False),
    'F': PaintingInfo(name="Modern/AI Landscape", artist="Unknown", nationality="N/A", year=None, is_french=False, before_1900=False),
}

result_string = ""
for label in sorted(paintings.keys()):
    painting = paintings[label]
    if painting.is_french and painting.before_1900:
        result_string += 'Y'
    else:
        result_string += 'N'

print(result_string)