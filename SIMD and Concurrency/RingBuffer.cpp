#include "stdafx.h"
#include "RingBuffer.h"


RingBuffer::RingBuffer(int size)
{
	ringSize=size;
	ringBuffer[0] = new int[size];
	ringBuffer[1] = new int[size];
	head=tail=0;
}


RingBuffer::~RingBuffer(void)
{
	delete[] ringBuffer[0];
	delete[] ringBuffer[1];
}


int RingBuffer::NextPosition(int pos) //returns the position after pos in the ring.
{
	return (pos+1)%ringSize; 
}

bool RingBuffer::IsRingEmpty()
{
	return head==tail; //returns whether the ring is empty
}

bool RingBuffer::IsRingFull() //returns whether the ring is full.
{
	return tail==NextPosition(head);
}

bool RingBuffer::AddToRing(int val[2])
{
	ringLock.lock();
	if (IsRingFull())
	{
		ringLock.unlock();
		return false;
	}
	(ringBuffer[0])[head] = val[0];
	(ringBuffer[1])[head] = val[1];
	head=NextPosition(head);
	ringLock.unlock();
	return true;
}

int *  RingBuffer::RemoveFromRing() 
{
	ringLock.lock();
	if (IsRingEmpty())
	{
		ringLock.unlock();
		return nullptr;
	}
	int send[2] = { ringBuffer[0][tail], ringBuffer[1][tail] };
	tail = NextPosition(tail);
	ringLock.unlock();
	return &send[0];
}


bool RingBuffer::IsRingEmptyCheck()   //Return whether ring is empty and no jobs are being processed. 
{
	ringLock.lock();
	bool res=IsRingEmpty();
	ringLock.unlock();
	return res;
}