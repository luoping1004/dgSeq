def ExcludeGene(small, Big):
    exclude = []
    for node in small:
        if node not in Big:
            exclude.append(node)
    for node in exclude:
        small.remove(node)
    return small
