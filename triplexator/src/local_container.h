#ifndef FBUSKE_APPS_TRIPLEXATOR_HEADER_LC_H
#define FBUSKE_APPS_TRIPLEXATOR_HEADER_LC_H

using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{    
	template <
	typename TFiber,
	typename TNeedle,
	typename TBitArray,
	typename TSeqNo,
	typename TOffset,
	typename TOptions
	>
	class LocalContainer {
	public:
		LocalContainer(
				TFiber &_alignedFiber,
				TNeedle &_alignedNeedle,
				TFiber const &_fiber,
				TNeedle const &_needle,
				TBitArray *_bitFiber,
				TSeqNo _fiberSeqNo,
				TSeqNo _needleSeqNo,
				TOffset _fiberPosOffset,
				unsigned const k,
				unsigned const alphabetSize,
				bool const plusStrand,
				TOptions const &_options);

		TFiber 		const	&fiber;
		TFiber 				&alignedFiber;
		TNeedle 	const	&needle;
		TNeedle 			&alignedNeedle;
		TBitArray			*bitFiber;
		TSeqNo		const	fiberSeqNo;
		TSeqNo		const	needleSeqNo;
		TOffset 			&fiberPosOffset;
		TOptions	const	&options;

		unsigned const k;
		unsigned const alphabetSize;
		bool const plusStrand;
	};

} //namespace SEQAN_NAMESPACE_MAIN

#endif
