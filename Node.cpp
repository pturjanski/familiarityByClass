#include "Node.h"

Node* Node::findChild(char c)
{
    for (unsigned int i=0;i<mChildren.size();i++){
        Node* tmp = mChildren.at(i);
        if (tmp->content()==c){
            return tmp;
        }
    }

    return NULL;
}
