#pragma once
#include <mutex>

class RingBuffer
{
private:
	int ringSize;
	int head, tail;
	int * ringBuffer[2];
	std::recursive_mutex ringLock;
private:
	int RingBuffer::NextPosition(int pos);
	bool IsRingFull();
	bool IsRingEmpty();
public:
	RingBuffer(int size);
	~RingBuffer(void);
	int * RemoveFromRing();
	bool AddToRing(int[2]);
	bool IsRingEmptyCheck();  //Return whether ring is empty
};
