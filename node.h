#ifndef NODE_H
#define NODE_H


class Node
{
    public:
        int index;
        double x;
        double y;


        Node();
        Node(int index, double x, double y);
        void draw(void);
};

#endif // NODE_H
