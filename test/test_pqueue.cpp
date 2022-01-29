#include <priority_queue.hpp>


int main(int argc, char ** argv)
{

	PriorityQueue<double> q(4);
	//add elements with priority
	q.Push(0,3.0);
	q.Push(1,5.0);
	q.Push(2,0.5);
	q.Push(3,10.0);
	//change priority
	q.DecreaseKey(3,2.0);
	q.IncreaseKey(1,9.0);
	//print contents (expect: 2 3 0 1)
	while(!q.Empty()) std::cout << q.Pop() << " ";
	std::cout << std::endl;
	return 0;
}