#ifndef _PRIORITY_QUEUE_H
#define _PRIORITY_QUEUE_H

#include <iostream>
#include <vector>
#include <string>
#include <limits.h>
#include <cmath>
#include <assert.h>


template<typename KeyType, typename Comparator = std::less<KeyType> >
class PriorityQueue
{
	int size_max, size_cur;
	std::vector<int> heap;
	std::vector<int> index;
	std::vector<KeyType> keys;
	void Swap(int i, int j)
	{
		int t = heap[i];
		heap[i] = heap[j];
		heap[j] = t;
		index[heap[i]] = i;
		index[heap[j]] = j;
	}
	void BubbleUp(int k)
	{
		while(k > 1 && !Comparator()(keys[heap[k/2]],keys[heap[k]]))
		{
			Swap(k, k/2);
			k = k/2;
		}
	}
	void BubbleDown(int k)
	{
		int j;
		while(2*k <= size_cur)
		{
			j = 2*k;
			if(j < size_cur && !Comparator()(keys[heap[j]],keys[heap[j+1]]) )
				j++;
			if(Comparator()(keys[heap[k]],keys[heap[j]]) )
				break;
			Swap(k, j);
			k = j;
		}
	}
public:
	PriorityQueue(int num_elements)
	{
		size_max = num_elements;
		size_cur = 0;
		keys.resize(size_max);
		heap.resize(size_max+1,INT_MAX);
		index.resize(size_max+1,INT_MAX);
	}
	void Push(int i, const KeyType & key)
	{
		size_cur++;
		index[i] = size_cur;
		heap[size_cur] = i;
		keys[i] = key;
		BubbleUp(size_cur);
	}
	int Pop()
	{
		if( size_cur == 0 ) 
			return INT_MAX;
		int min = heap[1];
		Swap(1, size_cur--);
		BubbleDown(1);
		index[min] = INT_MAX;
		heap[size_cur+1] = INT_MAX;
		return min;
	}
	int Peek() const
	{
		if( size_cur == 0 ) return INT_MAX;
		return heap[1];
	}
	void DecreaseKey(int i, const KeyType & key)
	{
		keys[i] = key;
		BubbleUp(index[i]);
	}
	void IncreaseKey(int i, const KeyType & key)
	{
		keys[i] = key;
		BubbleDown(index[i]);
	}
	void ChangeKey(int i, KeyType key)
	{
		keys[i] = key;
		if( !Comparator()(keys[i],key) )	
			BubbleUp(index[i]);
		else
			BubbleDown(index[i]);
	}
	const KeyType & GetKey(int i) const {return keys[i];}
	void Clear() {while( !Empty() ) Pop();}
	bool Contains(int i) const {return index[i] != INT_MAX;}
	int  Size() const {return size_cur;}
	bool Empty() const {return size_cur == 0;}
};

#endif
