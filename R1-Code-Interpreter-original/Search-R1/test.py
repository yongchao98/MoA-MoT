from typing import Union, List, Tuple, Dict, Optional
import re

def postprocess_predictions(predictions):
    """
    Process (text-based) predictions from llm into actions and validity flags.

    Args:
        predictions: List of raw predictions

    Returns:
        Tuple of (actions list, validity flags list)
    """
    actions = []
    contents = []

    for prediction in predictions:
        if isinstance(prediction, str): # for llm output
            pattern = r'<(search|answer)>(.*?)</\1>'
            match = re.search(pattern, prediction, re.DOTALL)
            if match:
                content = match.group(2).strip()  # Return only the content inside the tags
                action = match.group(1)
            else:
                content = ''
                action = None
        else:
            raise ValueError(f"Invalid prediction type: {type(prediction)}")

        actions.append(action)
        contents.append(content)

    return actions, contents

response_example_1 = '''
<think> I don't know where the first ICC World Cup was held. </think><search>when was the first icc world cup held</search>
<information>Doc 1(Title: "Cricket World Cup") the matches with coloured uniforms was the WSC Australians in wattle gold versus WSC West Indians in coral pink, played at VFL Park in Melbourne on 17 January 1979. The success and popularity of the domestic one-day competitions in England and other parts of the world, as well as the early One-Day Internationals, prompted the ICC to consider organising a Cricket World Cup. The inaugural Cricket World Cup was hosted in 1975 by England, the only nation able to put forward the resources to stage an event of such magnitude at the time. The 1975 tournament started on 7 June.

Doc 2(Title: "History of the ICC Cricket World Cup") first held in 1975 in England, which was the only nation that was able to put forward the resources to stage an event of such magnitude. The first three tournaments were officially known as the Prudential Cup after the financial services company Prudential plc who sponsored the event. The matches were of 60 overs per team and played in traditional white uniform and red balls. They were all played during the day. Eight teams participated in the first tournament: Australia, England, the West Indies, Pakistan, India, and New Zealand (the six Test nations at the time), with Sri Lanka and

'''

response_example_2 = '''
Doc 3(Title: "Cricket World Cup") Cricket World Cup The ICC Cricket World Cup is the international championship of One Day International (ODI) cricket. The event is organised by the sport's governing body, the International Cricket Council (ICC), every four years, with preliminary qualification rounds leading up to a finals tournament. The tournament is one of the world's most viewed sporting events and is considered the ""flagship event of the international cricket calendar"" by the ICC. The first World Cup was organised in England in June 1975, with the first ODI cricket match having been played only four years earlier. However, a separate Women's Cricket World</information>
<answer> 1979 </answer>
'''

actions, contents = postprocess_predictions([response_example_1, response_example_2])
print(actions)  # ['think', 'search', 'information', 'answer']
print(contents)