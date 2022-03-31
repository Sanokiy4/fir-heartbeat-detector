#pragma once
class RingBuffer
{
public:
	int size;
	RingBuffer();
private:
	float* buffer;
};

