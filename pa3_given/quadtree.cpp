/**
 * @file quadtree.cpp
 * Quadtree class implementation.
 */

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>



using namespace std;

#include "quadtree.h"
#include "png.h"



// Quadtree
//   - parameters: none
//   - constructor for the Quadtree class; makes an empty tree
Quadtree::Quadtree() 
{
    root = NULL;
    res = 0;
}

// Quadtree
//   - parameters: PNG const & source - reference to a const PNG
//                    object, from which the Quadtree will be built
//                 int resolution - resolution of the portion of source
//                    from which this tree will be built
//   - constructor for the Quadtree class; creates a Quadtree representing
//        the resolution by resolution block in the upper-left corner of
//        source
Quadtree::Quadtree(PNG const& source, int setresolution)
{
    
    root = NULL;
    this->buildTree(source, setresolution);
    

}

// Quadtree
//   - parameters: Quadtree const & other - reference to a const Quadtree
//                    object, which the current Quadtree will be a copy of
//   - copy constructor for the Quadtree class
Quadtree::Quadtree(Quadtree const& other)

{
    res = other.res;
    if(!other.root){
        root = NULL;
            }
    copy(other);
}


void Quadtree::copy(Quadtree const& other)
{
    res = other.res;
    
    root = new QuadtreeNode(other.root->element);
    root->level = log2(res);
    copy_helper(root , root->nwChild, other.root->nwChild);
    copy_helper(root , root->neChild, other.root->neChild);
    copy_helper(root , root->swChild, other.root->swChild);
    copy_helper(root , root->seChild, other.root->seChild);

    
}

void Quadtree::copy_helper(QuadtreeNode* &thisNodePrev,QuadtreeNode* &thisNode, QuadtreeNode* other)
{
    if (other == NULL)
        thisNode = NULL;
    

    else{
    thisNode = new QuadtreeNode(thisNodePrev ,other->element);
    copy_helper(thisNode,thisNode->swChild, other->swChild);
    copy_helper(thisNode,thisNode->seChild, other->seChild);
    copy_helper(thisNode,thisNode->nwChild, other->nwChild);
    copy_helper(thisNode,thisNode->neChild, other->neChild);
        
    }
    

    
}
// ~Quadtree
//   - parameters: none
//   - destructor for the Quadtree clas
Quadtree::~Quadtree()
{
    destroy_quadtree(root);
    res = 0;

}


void Quadtree::destroy_quadtree(QuadtreeNode *& node){



    if (node != NULL){

        
        destroy_quadtree(node->neChild);
        destroy_quadtree(node->nwChild);
        destroy_quadtree(node->seChild);
        destroy_quadtree(node->swChild);
    
    

        delete node;
        node = NULL;
    }
    

}

// operator=
//   - parameters: Quadtree const & other - reference to a const Quadtree
//                    object, which the current Quadtree will be a copy of
//   - return value: a const reference to the current Quadtree
//   - assignment operator for the Quadtree class
Quadtree const& Quadtree::operator=(Quadtree const& other)
{
    
    if(this != &other){
        this->~Quadtree();
        copy(other);
    }
    return *this;
        

}

// buildTree (public interface)
//   - parameters: PNG const & source - reference to a const PNG
//                    object, from which the Quadtree will be built
//                 int resolution - resolution of the portion of source
//                    from which this tree will be built
//   - transforms the current Quadtree into a Quadtree representing
//        the resolution by resolution block in the upper-left corner of
//        source
void Quadtree::buildTree(PNG const& source, int setresolution)
{

    this->~Quadtree();

    res = setresolution;

    root = new QuadtreeNode(this);
    buildTreeHelper(root, source, 0, 0, setresolution);



}



void Quadtree::buildTreeHelper(QuadtreeNode* node, PNG const& source, int x, int y, int res){

    

    if (node->level == 0){
        node->nwChild = NULL;
        node->neChild = NULL;
        node->swChild = NULL;
        node->seChild = NULL;
        node->element = *(source(x,y));
        return;
    }
    
    else {
        
        node->nwChild = new QuadtreeNode(node);
        node->neChild = new QuadtreeNode(node);
        node->swChild = new QuadtreeNode(node);
        node->seChild = new QuadtreeNode(node);
        
        buildTreeHelper(node->nwChild, source, x, y, res / 2);
        
        buildTreeHelper(node->neChild,  source, x + (res / 2), y, res / 2);
        
        buildTreeHelper(node->swChild, source, x, y + (res / 2), res / 2);
        
        buildTreeHelper(node->seChild, source, x + (res / 2), y + (res / 2), res / 2);
        
        
        
        
        node->element.red = ((node->nwChild)->element.red+(node->neChild)->element.red+(node->swChild)->element.red+(node->seChild)->element.red)/4;
        node->element.green = ((node->nwChild)->element.green+(node->neChild)->element.green+(node->swChild)->element.green+(node->seChild)->element.green)/4;
        node->element.blue = ((node->nwChild)->element.blue+(node->neChild)->element.blue+(node->swChild)->element.blue+(node->seChild)->element.blue)/4;
        node->element.alpha = ((node->nwChild)->element.alpha+(node->neChild)->element.alpha+(node->swChild)->element.alpha+(node->seChild)->element.alpha)/4;
        
    }

}

// getPixel (public interface)
//   - parameters: int x, int y - coordinates of the pixel to be retrieved
//   - return value: an RGBAPixel representing the desired pixel of the
//        underlying bitmap
//   - retrieves and returns the pixel at coordinates (x, y) in the
//        underlying bitmap
RGBAPixel Quadtree::getPixel(int x, int y) const
{
    if (x > res || y > res ) return RGBAPixel();
   
    if (!root) return RGBAPixel();
    
    else
        return getPixel_helper(root,  x,  y, res );
}


RGBAPixel Quadtree::getPixel_helper(QuadtreeNode* node,int  x,int y,int resolution) const{

    
    if(x< resolution/2 && y < resolution/2 && node->nwChild)
    {
        return getPixel_helper(node->nwChild, x, y, resolution/2);
    }
    
    else if(x>=resolution/2 && y >= resolution/2 && node->seChild)
    {
        return getPixel_helper(node->seChild, x-(resolution/2), y-(resolution/2), resolution/2);
    }
    
    else if(x<resolution/2 && y >= resolution/2 && node->swChild)
    {
        return getPixel_helper(node->swChild, x, y-(resolution/2), resolution/2);
    }
    
    else if (x>=resolution/2 && y < resolution/2 && node->seChild)
    {
        return getPixel_helper(node->neChild, x-(resolution/2), y, resolution/2);
    }
    
    else
        return node->element;
}


// decompress (public interface)
//   - parameters: none
//   - return value: a PNG object representing this quadtree's underlying
//        bitmap
//   - constructs and returns this quadtree's underlying bitmap
PNG Quadtree::decompress() const
{
    if(root){
    
    
    PNG source = PNG(res, res);
    
    
     decompress_helper(source);
        return source;}
    
    else
      return PNG();

}


void Quadtree::decompress_helper( PNG &source) const
{
    

            for(int i = 0 ; i<res; i++){
            for(int j = 0; j<res; j++){
                *(source(i,j)) = getPixel(i,j);
            }
        }
        
    
    

}


// clockwiseRotate (public interface)
//   - parameters: none
//   - transforms this quadtree into a quadtree representing the same
//        bitmap, rotated 90 degrees clockwise
void Quadtree::clockwiseRotate()
{

    clockwiseRotate_helper(root);
    
    
    
    
}


void Quadtree::clockwiseRotate_helper(QuadtreeNode * node){


    if (node->neChild){
        
        QuadtreeNode* tmpPtr;
    
        tmpPtr = node->neChild;
        node->neChild = node->nwChild;
        node->nwChild = node->swChild;
        node->swChild = node->seChild;
        node->seChild = tmpPtr;
        
        clockwiseRotate_helper(node->neChild);
        clockwiseRotate_helper(node->nwChild);
        clockwiseRotate_helper(node->seChild);
        clockwiseRotate_helper(node->swChild);
        
    }
    
    
}




// prune (public interface)
//   - parameters: int tolerance - an integer representing the maximum
//                    "distance" which we will permit between a node's color
//                    (i.e. the average of its descendant leaves' colors)
//                    and the color of each of that node's descendant leaves
//   - for each node in the quadtree, if the "distance" between the average
//        of that node's descendant leaves' colors and the color of each of
//        that node's descendant leaves is at most tolerance, this function
//        deletes the subtrees beneath that node; we will let the node's
//        color "stand in for" the colors of all (deleted) leaves beneath it
void Quadtree::prune(int tolerance)
{
    
    prune_helper(root, tolerance);
    
    
        }
    
    
        
        
        
void Quadtree::prune_helper(QuadtreeNode * node, int tolerance){
    

    if (node->nwChild){
    if (prunable (node ,node, tolerance)){

          destroy_quadtree(node->neChild);
          destroy_quadtree(node->nwChild);
          destroy_quadtree(node->swChild);
          destroy_quadtree(node->seChild);
        
        

        
    }
        
        
        else {
            
            prune_helper(node->nwChild, tolerance);
            prune_helper(node->neChild, tolerance);
            prune_helper(node->swChild, tolerance);
            prune_helper(node->seChild, tolerance);

        }
    }
    
    //cout << " we are done at prune_helper " << endl;
    
    
    

        
}
        
    

bool Quadtree::prunable(QuadtreeNode* root, QuadtreeNode* node, int tolerance) const{
    
    
    if(!node->nwChild){
    
        return (!(diff_of_Color(node->element, root->element) > tolerance));}
    
    
    else{

        
        return(
               prunable(root, node->nwChild, tolerance) &&
               
                prunable(root, node->neChild, tolerance) &&
               
               prunable(root, node->seChild, tolerance) &&
               
               prunable(root, node->swChild, tolerance));

    }
    
}



int Quadtree::diff_of_Color(RGBAPixel el, RGBAPixel el1) const{
    
    


    
    int diff = pow((el.red - el1.red),2) + pow((el.green - el1.green),2)+ pow((el.blue - el1.blue),2);
    

    
    return diff;
}

// pruneSize (public interface)
//   - parameters: int tolerance - an integer representing the maximum
//                    "distance" which we will permit between a node's color
//                    (i.e. the average of its descendant leaves' colors)
//                    and the color of each of that node's descendant leaves
//   - returns the number of leaves which this quadtree would contain if it
//        was pruned using the given tolerance; does not actually modify the
//        tree
int Quadtree::pruneSize(int tolerance) const
{
    int acm=0;
    //prunsize_calculate(root->level, acm1);
    //acm1++;
    
    int acm1 = pow(4, root->level);
    
    pruneSize_helper(tolerance, acm, root);

    return acm1-acm;
}


void Quadtree::pruneSize_helper(int tolerance, int & acm, QuadtreeNode * node) const{

    
    if (node->nwChild){
        if (prunable (node ,node, tolerance)){
            
            acm--;
            prunsize_leaves_calculate_acc(node->level, acm);
            
            
        }
        
        
        else {
            pruneSize_helper(tolerance, acm, node->nwChild);
            pruneSize_helper(tolerance, acm, node->neChild);
            pruneSize_helper(tolerance, acm, node->swChild);
            pruneSize_helper(tolerance, acm, node->seChild);
            
            
        }
    }

    


}



void Quadtree::prunsize_leaves_calculate_acc(int n, int & acm) const{
    
    acm += pow(4, n);

    
    
}




// idealPrune (public interface)
//   - parameters: int numLeaves - an integer representing the number of
//                    leaves we wish the quadtree to have, after pruning
//   - returns the minimum tolerance such that pruning with that tolerance
//        would yield a tree with at most numLeaves leaves
int Quadtree::idealPrune(int numLeaves) const
{
  


    return idealPrune_helper(numLeaves, 256*256/2, 256*256/4);

}





int Quadtree::idealPrune_helper(int numLeaves, int node, int a) const{

    
    if (pruneSize(node)> numLeaves){
        

        return idealPrune_helper(numLeaves, node + a, a/2);

    
    }
    
    if (pruneSize(node) < numLeaves)
    {
       

        return idealPrune_helper(numLeaves, node - a, a/2);
    
    }
    
    else {
        int b = node - 2*a;
        while(pruneSize(b) != numLeaves){
            b++;
        }
        
        return b;
    }

}

// QuadtreeNode
//   - parameters: none
//   - constructor for the QuadtreeNode class; creates an empty
//        QuadtreeNode, with all child pointers NULL
Quadtree::QuadtreeNode::QuadtreeNode()
{
    neChild = seChild = nwChild = swChild = NULL;
}
Quadtree::QuadtreeNode::QuadtreeNode(Quadtree * ob){

    neChild = seChild = nwChild = swChild = NULL;
    level = log2(ob->res);

}






Quadtree::QuadtreeNode::QuadtreeNode(QuadtreeNode const*  node){

    level = node->level -1;
    neChild = seChild = nwChild = swChild = NULL;
}


Quadtree::QuadtreeNode::QuadtreeNode(QuadtreeNode const*  node , RGBAPixel const& elem){
    element = elem;
    level = node->level -1;
    neChild = seChild = nwChild = swChild = NULL;
}




// QuadtreeNode
//   - parameters: RGBAPixel const & elem - reference to a const
//        RGBAPixel which we want to store in this node
//   - constructor for the QuadtreeNode class; creates a QuadtreeNode
//        with element elem and all child pointers NULL
Quadtree::QuadtreeNode::QuadtreeNode(RGBAPixel const& elem)
{
        element = elem;
    neChild = seChild = nwChild = swChild = NULL;
}


/**** Testing/grading functions ****/

// printTree (public interface)
//   - parameters: none
//   - prints the contents of the Quadtree using a preorder traversal
void Quadtree::printTree(ostream& out /* = cout */) const
{
    if (root == NULL)
        out << "Empty tree.\n";
    else{
        
        printTree(out, root, 1);
        
    }
}

// printTree (private helper)
//   - parameters: QuadtreeNode *current - pointer to the root of the
//                    subQuadtree which we wish to print
//                 int level - the current recursion depth; used for
//                    determining when to terminate recursion (see note below)
//   - prints the leaves of the Quadtree using a preorder traversal
void Quadtree::printTree(ostream& out, QuadtreeNode const* current,
                         int level) const
{
    // Is this a leaf?
    // Note: it suffices to check only one of the child pointers,
    // since each node should have exactly zero or four children.
    if(current->neChild == NULL) {
        out << current->element << " at depth " << level << "\n";
        return;
    }
    
    
    // This clause added for the sake of grading; we will never call
    // printTree() on quadtrees larger than  128x128.  (This is a
    // necessary restriction, as the grading scripts choke on programs
    // which produce excessive amounts of output.)
    if (level > 7) {
        out << "Error: infinite loop detected in printTree();"
        << " quadtree has a loop.\n";
        out << "Aborting program.\n";
        exit(1);
    }
    
    // Standard preorder traversal
    printTree(out, current->neChild, level + 1);
    printTree(out, current->seChild, level + 1);
    printTree(out, current->swChild, level + 1);
    printTree(out, current->nwChild, level + 1);
}

// operator==
//   - parameters: Quadtree const & other - reference to a const Quadtree
//                    object, against which the current Quadtree will be
//                    compared
//   - return value: a boolean which is true if the Quadtrees are deemed
//        "equal", and false otherwise
//   - compars the current Quadtree with the parameter Quadtree, and
//        determines whether or not the two are the same
// Note: this method relies on the private helper method compareTrees()
bool Quadtree::operator==(Quadtree const& other) const
{
    return compareTrees(root, other.root);
}

// compareTrees
//   - parameters: QuadtreeNode const * firstPtr - pointer to the root
//                    of a subtree of the "first" Quadtree under
//                    consideration
//                 QuadtreeNode const * secondPtr - pointer to the root
//                    of a subtree of the "second" Quadtree under
//                    consideration
//   - return value: a boolean which is true if the subQuadtrees are deemed
//        "equal", and false otherwise
//   - compares the subQuadtree rooted at firstPtr with the subQuadtree
//        rooted at secondPtr, and determines whether the two are the same
//   - this function only compares the leaves of the trees, as we did not
//     impose any requirements on what you should do with interior nodes
bool Quadtree::compareTrees(QuadtreeNode const* firstPtr,
                            QuadtreeNode const* secondPtr) const
{
    if (firstPtr == NULL && secondPtr == NULL)
        return true;
    
    if (firstPtr == NULL || secondPtr == NULL)
        return false;
    
    // if they're both leaves, see if their elements are equal
    // note: child pointers should _all_ either be NULL or non-NULL,
    // so it suffices to check only one of each
    if (firstPtr->neChild == NULL && secondPtr->neChild == NULL) {
        if (firstPtr->element.red != secondPtr->element.red
            || firstPtr->element.green != secondPtr->element.green
            || firstPtr->element.blue != secondPtr->element.blue)
            return false;
        else
            return true;
    }
    
    // they aren't both leaves, so recurse
    return (compareTrees(firstPtr->neChild, secondPtr->neChild)
            && compareTrees(firstPtr->nwChild, secondPtr->nwChild)
            && compareTrees(firstPtr->seChild, secondPtr->seChild)
            && compareTrees(firstPtr->swChild, secondPtr->swChild));
}








