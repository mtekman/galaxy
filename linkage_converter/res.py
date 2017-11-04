

marker_order = [ "rs123", "rs111111119", "rs15441", "rs1" ]

def __generateHeaders(marker_order, npad_left = 10):
    max_len = -1
    markerpadd = []
    for marker in marker_order:
        nmark = len(marker)
        if nmark > max_len:
            max_len = nmark

    # paddleft
    for marker in marker_order:
        markerpadd.append(("%%-%ds" % max_len) % marker)

    # transpose
    buffer_left = ("%%%ds" % npad_left) % " "    
    return '\n'.join([buffer_left + "  ".join(x) for x in zip(*markerpadd)][::-1])


print(__generateHeaders(marker_order, 30))
