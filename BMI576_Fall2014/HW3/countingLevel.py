# """
# Counting level of a node in a tree
# """


#    import pdb; pdb.set_trace()

#inv_relations = {1: [2, 5], 2: [3, 4]}
#root = 1

class BinTree:
  """Node in a binary tree
##       1
##     2     3
##   4   5  6  7
##  8
  treeRep = [(1,2,3),(2,4,5),(3,6,7),(4,8,None)]
  tree= BinTree.createTree(treeRep)
  tree.printBfsLevels()
>>>
1 

2 3 

4 5 6 7 

"""
  def __init__(self,val,leftChild=None,rightChild=None,root=None):
    self.val=val
    self.leftChild=leftChild
    self.rightChild=rightChild
    self.root=root
    if not leftChild and not rightChild:
      self.isExternal=True

  def getChildren(self,node):
    children=[]
    if node.isExternal:
      return []
    if node.leftChild:
      children.append(node.leftChild)
    if node.rightChild:
      children.append(node.rightChild)
    return children

  @staticmethod
  def createTree(tupleList):
    "Creates a Binary tree Object from a given Tuple List"
    Nodes={}
    root=None
    for item in treeRep:
      if not root:
        root=BinTree(item[0])
        root.isExternal=False
        Nodes[item[0]]=root
        root.root=root
        root.leftChild=BinTree(item[1],root=root)
        Nodes[item[1]]=root.leftChild
        root.rightChild=BinTree(item[2],root=root)
        Nodes[item[2]]=root.rightChild
      else:
        CurrentParent=Nodes[item[0]]
        CurrentParent.isExternal=False
        CurrentParent.leftChild=BinTree(item[1],root=root)
        Nodes[item[1]]=CurrentParent.leftChild
        CurrentParent.rightChild=BinTree(item[2],root=root)
        Nodes[item[2]]=CurrentParent.rightChild
    root.nodeDict=Nodes
    return root

  def printBfsLevels(self, levels=None):
#      import pdb; pdb.set_trace()
      global store
      global index
      if levels==None:
          levels=[self]
      nextLevel=[]
      for node in levels:
          store.insert(index, node.val)
          index += 1
          print node.val,
      for node in levels:
          nextLevel.extend(node.getChildren(node))
      print '\n'
      if nextLevel:
          node.printBfsLevels(nextLevel)
    
if __name__ == '__main__':
    store = []
    index = 0
    level = []
    treeRep = [(1,2,5),(2,3,4)]
    tree = BinTree.createTree(treeRep)
    tree.printBfsLevels()
