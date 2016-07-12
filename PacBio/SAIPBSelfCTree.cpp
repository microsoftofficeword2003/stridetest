///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#include "SAIPBSelfCTree.h"

#include "BWTAlgorithms.h"
#include "stdaln.h"

using namespace std;


//
// Class: SAIPBSelfCorrectTree
SAIPBSelfCorrectTree::SAIPBSelfCorrectTree(
					const BWT* pBWT,
                    const BWT* pRBWT,
					size_t min_SA_threshold,
					int maxLeavesAllowed) :
                    m_pBWT(pBWT), m_pRBWT(pRBWT), m_min_SA_threshold(min_SA_threshold), m_maxLeavesAllowed(maxLeavesAllowed),
					m_expectedLength(0), m_currentLength(0), m_pRootNode(NULL), m_isSourceRepeat(false), 
					m_isTargetRepeat(false), m_isLargeLeaveRemoved(false), m_prevExtPos(0)
{
	// initialize google dense hashmap
	kmerHash.set_empty_key("");
	kmerHash.resize(6000);
}

//
SAIPBSelfCorrectTree::~SAIPBSelfCorrectTree()
{
    // Recursively destroy the tree
	if(m_pRootNode!=NULL)
		delete m_pRootNode;
		
	kmerHashiter hashiter = kmerHash.begin();
	for(; hashiter != kmerHash.end(); ++hashiter)
	{
		if(hashiter->second!=NULL)
		{
			delete hashiter->second;
			hashiter->second = NULL;
		}
	}

}

// create root node using src
void SAIPBSelfCorrectTree::initializeSearchTree(std::string src, size_t hashKmerSize)
{
	m_leaves.clear();

	// Create the root node containing the seed string
    SAIntervalNode *pRootNode = new SAIntervalNode(&src, NULL);	
	// store initial str of root
    pRootNode->computeInitial(src);
	
	// extract suffix of length minOverlap
    std::string beginningkmer=src.substr(src.length() - hashKmerSize);
	// std::cout << beginningkmer << "\n";

	// Compute initial SA intervals in both forward and reverse BWT
    pRootNode->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
    pRootNode->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(beginningkmer));
	m_leaves.push_back(pRootNode);
	
	// for cleaning memory
	m_pRootNode = pRootNode;	
	m_seedLength = src.length();
	m_currentLength = src.length();
	
	// debug
	// std::cout << beginningkmer << "\t" << pRootNode->fwdInterval.size() << "\t" << pRootNode->rvcInterval.size() << "\n";
	
}

void SAIPBSelfCorrectTree::initializeTerminalIntervals(std::string dest, size_t hashKmerSize)
{
	//ending kmer is a prefix of dest
    //initialize the ending SA intervals with kmer length=m_minOverlap	
	std::string endingkmer = dest.substr(0, hashKmerSize);
	m_fwdTerminatedInterval = BWTAlgorithms::findInterval( m_pRBWT, reverse(endingkmer));
	m_rvcTerminatedInterval = BWTAlgorithms::findInterval( m_pBWT, reverseComplement(endingkmer));
}

// find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
int SAIPBSelfCorrectTree::mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
											size_t hashKmerSize, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength)
{	
	initializeSearchTree(src, hashKmerSize);
	initializeTerminalIntervals(dest, hashKmerSize);
	m_expectedLength = expectedLength;
	size_t maxUsedLeaves=0;
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		refineSAInterval(hashKmerSize-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize,1);

		if(m_leaves.size()> maxUsedLeaves) maxUsedLeaves=m_leaves.size();
		
		// extension succeed
		if(!newLeaves.empty())
			m_currentLength++;

		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if( (size_t) m_currentLength >= minLength)
			isTerminated(results);

		// if(maxLength==100)
			printLeaves(hashKmerSize);
	}
	
	return findBestExtSeq(results, src, dest, mergedseq, hashKmerSize, expectedLength, maxLength, maxLeaves);
}

int SAIPBSelfCorrectTree::findBestExtSeq(SAIntervalNodeResultVector& results, const std::string &src, const std::string &dest, std::string &mergedseq,
										size_t hashKmerSize, size_t expectedLength, size_t maxLength, size_t maxLeaves)
{
	if(results.size()>0)
	{
		double maxKmerCoverage = 0;
		int minLengthDiff = 100000;

		for (size_t i = 0 ; i < results.size() ;i++)
		{
			std::string tmpseq;
			// bug fix: dest may be shorter than hashKmerSize
			if(dest.length() > hashKmerSize)
				tmpseq=results[i].thread + dest.substr(hashKmerSize);
			else
				tmpseq=results[i].thread;

			int currLengthDiff = std::abs((int)tmpseq.length()-(int)expectedLength);

			/**** Math bug: the division requires denominator larger than numerator for accurate average frequency
					numerator denominator division
					100	1000	0.1
					99	999	0.099099099
					100	10	10
					99	9	11
			***/			
			double avgCov = (double)results[i].SAICoverage /(tmpseq.length()+1000000);

			// std::cout << ">" << currLengthDiff << "\t" << avgCov <<  "\t" << tmpseq.length()-src.length()-dest.length() << "\t" << expectedLength << "\n";
			
			// if(maxLength==97 || maxLength==152)
				// std::cout << ">" << src.length() << "\n" << tmpseq.substr(src.length()) <<  "\n";
			
			bool isLengthDiffBetter = currLengthDiff < minLengthDiff && std::abs(currLengthDiff - minLengthDiff) >3;
			
			// But the length diff is not informative for small indels, where kmer coverage is more informative			
			bool isKmerCoverageBetter = std::abs(currLengthDiff - minLengthDiff) <=3 && (maxKmerCoverage<avgCov);
			if(isLengthDiffBetter || isKmerCoverageBetter)
			{
				minLengthDiff = currLengthDiff;
				maxKmerCoverage=avgCov;
				mergedseq=tmpseq;
			}
		}
		
		// std::cout << ">Selected: " << results.size() << "\t" << maxKmerCoverage << "\t" << minLengthDiff  <<"\n";
		// std::cout << mergedseq.substr(src.length()-hashKmerSize) << "\n";
		// std::cout << mergedseq.length()-src.length() << "\n";

		return 1;
	}
				
    // Did not reach the terminal kmer
    if(m_leaves.empty() && m_currentLength >= (int)(expectedLength-m_seedLength)/2+m_seedLength)
        return -1;	// high error
    else if( (size_t) m_currentLength > maxLength)
        return -2;	// exceed search depth
    else if(m_leaves.size() > maxLeaves)
        return -3;	// too much repeats
	else if(m_leaves.empty() && m_currentLength < (int)(expectedLength-m_seedLength)/2+m_seedLength)
		return -4;
	else
		return -5;
}
// find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
int SAIPBSelfCorrectTree::extOverUnique(const std::string &src, const std::string &dest, std::string &mergedseq, 
											size_t hashKmerSize, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength)
{	
	initializeSearchTree(src, hashKmerSize);
	initializeTerminalIntervals(dest, hashKmerSize);
	m_expectedLength = expectedLength;

	size_t maxUsedLeaves=0;
	
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		refineSAInterval(hashKmerSize-1);
		
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize,1);

		if(m_leaves.size()> maxUsedLeaves) maxUsedLeaves=m_leaves.size();
		
		// extension succeed
		if(newLeaves.empty() || m_isLargeLeaveRemoved)
		{
			// current nodes were changed
			if(m_isLargeLeaveRemoved && !newLeaves.empty())
			{
				rollBackExtension(hashKmerSize-1);
				newLeaves.clear();
			}

			extendByNewSeeds(src,  dest, hashKmerSize, maxLength);	// add new seeds from extended sequences			
			attempToExtendUsingHash(newLeaves, hashKmerSize,1);
			m_isLargeLeaveRemoved = false;	// reset flag
		}
		
		if(!newLeaves.empty())
			m_currentLength++;
		
		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if( (size_t) m_currentLength >= minLength)
			isTerminated(results);

		// if(maxLength==190)
			printLeaves(hashKmerSize);
		
	}
	
	return findBestExtSeq(results, src, dest, mergedseq, hashKmerSize, expectedLength, maxLength, maxLeaves);

}

// find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
int SAIPBSelfCorrectTree::extOverRepeat(const std::string &src, const std::string &dest, std::string &mergedseq, 
											size_t hashKmerSize, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength,std::string currseq)
{	
	size_t globalKmerSize = 18;
	initializeSearchTree(src, hashKmerSize);
	initializeTerminalIntervals(dest, hashKmerSize);
	m_expectedLength = expectedLength;
	m_min_SA_threshold = 2;
	size_t maxUsedLeaves=0;
	SAIntervalNodeResultVector results;
	while(!m_leaves.empty() && m_leaves.size() <= maxLeaves &&  (size_t)m_currentLength <= maxLength)
	{
		// std::cout << m_leaves.size() << "\t" << currentLength << "\t" << expectedLength << "\t" << maxLength << "\n";
		if(m_currentLength<18)
			globalKmerSize = hashKmerSize;
		else
			globalKmerSize = 18;
			
		refineSAInterval(globalKmerSize-1);
			 
		STNodePtrList newLeaves;
	
		// attempt to extend one base for each leave
		// ACGT-extend the leaf nodes via updating existing SA interval
		attempToExtendUsingHash(newLeaves, hashKmerSize, 2); 

		if(m_leaves.size()> maxUsedLeaves) maxUsedLeaves=m_leaves.size();
		
		// extension succeed
		if(newLeaves.empty() || m_isLargeLeaveRemoved)
		{			
			// current nodes were changed
			if(m_isLargeLeaveRemoved && !newLeaves.empty())
			{
				rollBackExtension(globalKmerSize-1);
				newLeaves.clear();
			}
			cout <<m_leaves.size()<<"leavesize\n";
			m_isLargeLeaveRemoved = false;	// reset flag
			NewextendByNewSeeds(src,  dest, hashKmerSize, maxLength,currseq);	// add new seeds from extended sequences
			refineSAInterval(globalKmerSize-1);
			attempToExtendUsingHash(newLeaves, hashKmerSize, 2);
			std::cout << "Level 1\n";
			
			// reduce global kmer freq threshold
			if(newLeaves.empty())
			{
				refineSAInterval(globalKmerSize-1);
				attempToExtendUsingHash(newLeaves, hashKmerSize, 1);
				std::cout << "Level 2\n";
				
				// reduce local kmer freq threshold
				if(newLeaves.empty())
				{
					m_min_SA_threshold = 1;
					refineSAInterval(globalKmerSize-1);
					attempToExtendUsingHash(newLeaves, hashKmerSize, 1);
					m_min_SA_threshold = 2;
					std::cout << "Level 3\n";
				}
			}

		}
		
		if(!newLeaves.empty())
			m_currentLength++;
		
		m_leaves.clear();
		m_leaves = newLeaves;
		
		// see if terminating string is reached
		if( (size_t) m_currentLength >= minLength)
			isTerminated(results);

		// if(maxLength==100)
			printLeaves(hashKmerSize);
		
	}
	
	return findBestExtSeq(results, src, dest, mergedseq, hashKmerSize, expectedLength, maxLength, maxLeaves);

}

void SAIPBSelfCorrectTree::NewextendByNewSeeds(const std::string& src, const std::string& dest, size_t hashKmerSize, size_t maxLength,std::string currseq)
{
	// search for the leave with highest avg freq
	double maxAvgCov = 0;
	SAIntervalNode* pMaxNode = m_leaves.front();
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		double currAvgCov = (double)(*iter)->getKmerCount() / (m_currentLength);
		// std::cout << currAvgCov << "\n";
		if( currAvgCov > maxAvgCov)
		{
			pMaxNode = *iter;
			maxAvgCov = currAvgCov;
		}
		else if (currAvgCov < maxAvgCov)
		{
			m_leaves.erase(iter);
			iter--;
		}
	}

	// std::cout << m_leaves.size() << "\t" << m_currentLength << "\t" << m_prevExtPos <<" add new kmers\n";
    
	
	std::string extendedSeq = pMaxNode->getFullString();
	int matchLen = 0; 
	AlnAln *aln_global;
	aln_global = aln_stdaln(extendedSeq.c_str() ,currseq.substr(0,extendedSeq.length()).c_str(), &aln_param_blast, 1, 1);
	for(int i = 0 ; aln_global->outm[i] != '\0' ; i++)
				if(aln_global->outm[i] == '|')
					matchLen++;
	double matchPercent = (double)matchLen /extendedSeq.length();
	std::cout<< currseq.substr(0,extendedSeq.length()) << "|" << extendedSeq << "||" <<matchLen<<"+"<<extendedSeq.length()<<" gogopower" << endl;
	
	// size_t largeKmerSize = std::min(src.length(), dest.length());
	// largeKmerSize should be large enough for reducing false seeds
	size_t largeKmerSize = 17;

	// counter of consecutive errors
	// static size_t errCount = 0;

	if(extendedSeq.length() > largeKmerSize && extendedSeq.length() > src.length())
	{
		// errCount++;
		
		// the new seeds should not contribute too much new kmers, we limit to half of largeKmerSize
		addHashByExtSeed(extendedSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*0.5, maxLength);
		// addHashByExtSeed(extendedSeq, src.length(), largeKmerSize, hashKmerSize, maxLength);
		
		// don't repetitively collect new seeds for error extension
		if((int)m_currentLength - (int)m_prevExtPos >= 3 && (int)m_currentLength - (int)m_prevExtPos <= 9 )
		{	
			// std::cout << "haha " << m_currentLength << "\t" << m_prevExtPos << "\n";
			if((int)m_currentLength - (int)m_prevExtPos > 3)
				m_prevExtPos++;

			return;
		}
		
		if((int)m_currentLength - (int)m_prevExtPos > 9)
			m_prevExtPos = m_currentLength;

		// get second seed
		if(extendedSeq.length() > largeKmerSize*2.5)
		{
			std::string tmpSeq = extendedSeq.substr(0, extendedSeq.length()-largeKmerSize*0.5);
			addHashByExtSeed(tmpSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*1.5, maxLength);
		}
		
		// get third seed, the ext distance should be larger for better accuracy, the reason is probably due to indel errors which disrupt expected distance
		if(extendedSeq.length() > largeKmerSize*3.5)
		{
			std::string tmpSeq = extendedSeq.substr(0, extendedSeq.length()-largeKmerSize*1.5);
			addHashByExtSeed(tmpSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*2.5, maxLength);
		}
	}
}


void SAIPBSelfCorrectTree::extendByNewSeeds(const std::string& src, const std::string& dest, size_t hashKmerSize, size_t maxLength)
{
	// search for the leave with highest avg freq
	double maxAvgCov = 0;
	SAIntervalNode* pMaxNode = m_leaves.front();
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		double currAvgCov = (double)(*iter)->getKmerCount() / (m_currentLength);
		// std::cout << currAvgCov << "\n";
		if( currAvgCov > maxAvgCov)
		{
			pMaxNode = *iter;
			maxAvgCov = currAvgCov;
		}
		else if (currAvgCov < maxAvgCov)
		{
			m_leaves.erase(iter);
			iter--;
		}
	}

	// std::cout << m_leaves.size() << "\t" << m_currentLength << "\t" << m_prevExtPos <<" add new kmers\n";
    
	
	std::string extendedSeq = pMaxNode->getFullString();
	
	// size_t largeKmerSize = std::min(src.length(), dest.length());
	// largeKmerSize should be large enough for reducing false seeds
	size_t largeKmerSize = 17;

	// counter of consecutive errors
	// static size_t errCount = 0;

	if(extendedSeq.length() > largeKmerSize && extendedSeq.length() > src.length())
	{
		// errCount++;
		
		// the new seeds should not contribute too much new kmers, we limit to half of largeKmerSize
		addHashByExtSeed(extendedSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*0.5, maxLength);
		// addHashByExtSeed(extendedSeq, src.length(), largeKmerSize, hashKmerSize, maxLength);
		
		// don't repetitively collect new seeds for error extension
		if((int)m_currentLength - (int)m_prevExtPos >= 3 && (int)m_currentLength - (int)m_prevExtPos <= 9 )
		{	
			// std::cout << "haha " << m_currentLength << "\t" << m_prevExtPos << "\n";
			if((int)m_currentLength - (int)m_prevExtPos > 3)
				m_prevExtPos++;

			return;
		}
		
		if((int)m_currentLength - (int)m_prevExtPos > 9)
			m_prevExtPos = m_currentLength;

		// get second seed
		if(extendedSeq.length() > largeKmerSize*2.5)
		{
			std::string tmpSeq = extendedSeq.substr(0, extendedSeq.length()-largeKmerSize*0.5);
			addHashByExtSeed(tmpSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*1.5, maxLength);
		}
		
		// get third seed, the ext distance should be larger for better accuracy, the reason is probably due to indel errors which disrupt expected distance
		if(extendedSeq.length() > largeKmerSize*3.5)
		{
			std::string tmpSeq = extendedSeq.substr(0, extendedSeq.length()-largeKmerSize*1.5);
			addHashByExtSeed(tmpSeq, src.length(), largeKmerSize, hashKmerSize, extendedSeq.length()+largeKmerSize*2.5, maxLength);
		}
	}
}


// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
// Contaminated reads often are simple repeats C* or T* with large freq
// Give up this read if the freq is way too large
size_t SAIPBSelfCorrectTree::addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, bool skipRepeat, int expectedLength)
{
	// PacBio errors create repeat-like seeds with large interval
	// limit the upper bound of interval for speedup
	const int64_t maxIntervalSize = 40;
	
	// LF-mapping of each fwd index using largeKmerSize
	// assert(seedStr.length() >= largeKmerSize);
	std::string initKmer = seedStr.substr(seedStr.length() - largeKmerSize);
	BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));
    
	size_t kmerFreq = 0;
	kmerFreq += fwdInterval.isValid()?fwdInterval.size():0;
	kmerFreq += fwdInterval.isValid()?rvcInterval.size():0;

	// std::cout << initKmer << "\t" << smallKmerSize << "\t" << reverseComplement(initKmer) << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() << "\t" << kmerFreq << "\n";

	// skip repeat only in the 1st round
	if(skipRepeat && kmerFreq > 128) return kmerFreq;
	
	// extend each SA index and collect kmers of smallKmerSize along the extension
	for(int64_t fwdRootIndex = fwdInterval.lower; 
		fwdInterval.isValid() && fwdRootIndex <= fwdInterval.upper &&  fwdRootIndex - fwdInterval.lower < maxIntervalSize; 
		fwdRootIndex++)
	{
		// extract small hash Kmer
		std::string currentFwdKmer = seedStr.substr(seedStr.length() - smallKmerSize);

		// Bug fix: the first and last kmers in the seeds must be added.
		insertKmerToHash(currentFwdKmer, seedStr.length(), seedStr.length(), smallKmerSize, maxLength, expectedLength);		
		
		// extend the fwdIndex via LF mapping
		int64_t fwdIndex = fwdRootIndex;
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;
			
			// if(seedStr.length()==105)
				// std::cout << currentLength << ":" << currentFwdKmer << "\n";

			insertKmerToHash(currentFwdKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	
	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; 
		rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && rvcRootIndex - rvcInterval.lower < maxIntervalSize; 
		rvcRootIndex++)
	{
		std::string currentRvcKmer = reverseComplement( seedStr.substr(seedStr.length() - smallKmerSize) );
		insertKmerToHash(currentRvcKmer, seedStr.length(), seedStr.length(), smallKmerSize, maxLength, expectedLength);
		
		int64_t rvcIndex = rvcRootIndex;
		
		for(int64_t currentLength = (int64_t)seedStr.length(); currentLength <= (int64_t)maxLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, smallKmerSize-1);
			insertKmerToHash(currentRvcKmer, seedStr.length(), currentLength, smallKmerSize, maxLength, expectedLength);
			// if(seedStr.length()==105)
				// std::cout << currentLength << ":" << currentRvcKmer << "\n";

			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	// if(seedStr.length()==105)
	// {
		// kmerHashiter iter = kmerHash.find("GGCAAACCATG");/*GATGCGACGCTGTCG*/
		// if(iter != kmerHash.end()){
			// iter->second->printFreq();
			// // getchar();
		// }
		// std::cout << "Reverse Complement\n";
		// iter = kmerHash.find("CATGGTTTGCC");/*CGACAGCGTCGCATC*/
		// if(iter != kmerHash.end()){
			// iter->second->printFreq();
			// // getchar();
		// }
		
		// std::cout << iter1->second->getSumOfFreq(seedStr.length() - m_seedLength);
		// std::cout << kmerHash.size() << "\n";
	// }
	return kmerFreq;
}


size_t SAIPBSelfCorrectTree::addHashByExtSeed(std::string& extStr, size_t srcLength, size_t largeKmerSize, size_t smallKmerSize, size_t extLength, size_t maxLength)
{
	// PacBio errors create repeat-like seeds with large interval
	// limit the upper bound of interval for speedup
	const int64_t maxIntervalSize = 40;
	
	// LF-mapping of each fwd index using largeKmerSize
	
	size_t currkmersize=largeKmerSize;
	size_t kmerFreq = 0;
	std::string initKmer = extStr.substr(extStr.length() - currkmersize );
	BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));
	
	while (1)
	{
		initKmer = extStr.substr(extStr.length() - currkmersize );
			
		 fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(initKmer));
		 rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(initKmer));
		
		kmerFreq = 0;
		kmerFreq += fwdInterval.isValid()?fwdInterval.size():0;
		kmerFreq += fwdInterval.isValid()?rvcInterval.size():0; 
		cout<<initKmer<<"        ||      "<<kmerFreq<<"freqqqqqqqqqqqqqqqqqqqqqq\n";
		// std::cout << initKmer << "\t" << smallKmerSize << "\t" << reverseComplement(initKmer) << "\t" << fwdInterval.size() << "\t" << rvcInterval.size() << "\t" << kmerFreq << "\n";
		// if(isLowComplexity(initKmer, 0.9, true)) return 0;

		// skip repeat only in the 1st round
		if(kmerFreq > 64 && extStr.length()>currkmersize) currkmersize++;
		else if (kmerFreq == 0 && currkmersize > 13) currkmersize--;
		else if(kmerFreq > 64 ) return kmerFreq;
		else if (kmerFreq == 0 && currkmersize == 13) return kmerFreq;
		else break;
	}
	// if(kmerFreq > 64 ) return kmerFreq;
	
	// extend each SA index and collect kmers of smallKmerSize along the extension
	for(int64_t fwdRootIndex = fwdInterval.lower; 
		fwdInterval.isValid() && fwdRootIndex <= fwdInterval.upper &&  fwdRootIndex - fwdInterval.lower < maxIntervalSize; 
		fwdRootIndex++)
	{
		// extract small hash Kmer
		std::string currentFwdKmer = extStr.substr(extStr.length() - smallKmerSize);

		// Bug fix: the first and last kmers in the seeds must be added.
		insertKmerToHash(currentFwdKmer, srcLength, extStr.length(), smallKmerSize, maxLength, -1);
		
		// extend the fwdIndex via LF mapping
		int64_t fwdIndex = fwdRootIndex;
		for(int64_t currentLength = (int64_t)extStr.length(); currentLength <= (int64_t)extLength; currentLength++)
		{
			char b = m_pRBWT->getChar(fwdIndex);
			if(b == '$') break;
			
			currentFwdKmer = currentFwdKmer.substr(1) + b;
			
			// if(extStr.length()==105)
				// std::cout << currentLength << ":" << currentFwdKmer << "\n";

			insertKmerToHash(currentFwdKmer, srcLength, currentLength+1, smallKmerSize, maxLength, -1);

			// LF mapping
            fwdIndex = m_pRBWT->getPC(b) + m_pRBWT->getOcc(b, fwdIndex - 1);			
		}
	}

	
	// LF-mapping of each rvc index	
	for(int64_t rvcRootIndex=rvcInterval.lower; 
		rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid() && rvcRootIndex - rvcInterval.lower < maxIntervalSize; 
		rvcRootIndex++)
	{
		std::string currentRvcKmer = reverseComplement( extStr.substr(extStr.length() - smallKmerSize) );
		insertKmerToHash(currentRvcKmer, srcLength, extStr.length(), smallKmerSize, maxLength, -1);
		
		int64_t rvcIndex = rvcRootIndex;
		
		for(int64_t currentLength = (int64_t)extStr.length(); currentLength <= (int64_t)extLength; currentLength++)
		{
			char b = m_pBWT->getChar(rvcIndex);
			if(b == '$') break;
			
			currentRvcKmer = b + currentRvcKmer.substr(0, smallKmerSize-1);
			insertKmerToHash(currentRvcKmer, srcLength, currentLength+1, smallKmerSize, maxLength, -1);
			// if(extStr.length()==105)
				// std::cout << currentLength << ":" << currentRvcKmer << "\n";

			// LF mapping
            rvcIndex = m_pBWT->getPC(b) + m_pBWT->getOcc(b, rvcIndex - 1);
		}
	}
	
	// if(extStr.length()==105)
	// {
		// kmerHashiter iter = kmerHash.find("CAGCGTGGTGT");/*GATGCGACGCTGTCG*/
		// if(iter != kmerHash.end()){
			// iter->second->printFreq();
		// }
		// std::cout << "Reverse Complement\n";
		// iter = kmerHash.find("ACACCACGCTG");/*CGACAGCGTCGCATC*/
		// if(iter != kmerHash.end()){
			// iter->second->printFreq();
		// }

		// std::cout << kmerHash.size() << "\n";
	// }
	return kmerFreq;
}

void SAIPBSelfCorrectTree::insertKmerToHash(std::string& insertedKmer, size_t seedStrLen, size_t currentLength, size_t smallKmerSize, size_t maxLength, int expectedLength)
{
	kmerHashiter hashiter = kmerHash.find(insertedKmer);
	
	if(hashiter == kmerHash.end())
	{
		KmerFeatures *newEntry;
		// source to target
		if(expectedLength<0)
			newEntry = new KmerFeatures(currentLength - seedStrLen, maxLength);
		// target to source
		else
			newEntry = new KmerFeatures(expectedLength - currentLength + smallKmerSize, maxLength);
			
		kmerHash.insert(std::make_pair<std::string, KmerFeatures*>(insertedKmer, newEntry));
	}
	else
	{
		if(expectedLength<0)
			hashiter->second->add(currentLength - seedStrLen);
		else
			hashiter->second->add(expectedLength - currentLength + smallKmerSize);
	}
}

void SAIPBSelfCorrectTree::printLeaves(size_t hashKmerSize)
{
	std::cout << m_leaves.size() << ":" << m_currentLength << "\n";
	// cout << "CATTTATCTGAGATATAAAAAAACCCGCCCCTGACAGGCGGGAAGAACGGCAACTAAACTGTTATTCAGTGGCATTTAGATCTATGACGTATCTGGCAAACCATGCCCGATGCGACGCTGTCGCGT\n";
	for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		std::string STNodeStr = (*iter)->getFullString();
		std::string fwdKmer = (*iter)->getSuffix(hashKmerSize);
		kmerHashiter hashiter = kmerHash.find(fwdKmer);
		std::cout << STNodeStr.substr(m_seedLength-hashKmerSize);
		
		if(hashiter!=kmerHash.end() 
			// && kmerHash[fwdKmer]->getTotalFreq() > 0 
			// && kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0 
		  )
		{
			std::cout << " " << fwdKmer 
			// << ":" << kmerHash[fwdKmer]->getTotalFreq()
			// local kmer frequency at m_currentLength-m_seedLength
			<< ":" << kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength);
			// << ":" << kmerHash[fwdKmer]->getTotalSum()/kmerHash[fwdKmer]->getTotalFreq() + src.length()

			if(kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0)
				std::cout << ":" << kmerHash[fwdKmer]->getSumOfPos(m_currentLength-m_seedLength)/kmerHash[fwdKmer]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength;
		}
		
		std::string rvcKmer = reverseComplement(fwdKmer);
		kmerHashiter hashiter2 = kmerHash.find(rvcKmer);
		
		
		if(hashiter2!=kmerHash.end() 
			// && kmerHash[rvcKmer]->getTotalFreq() > 0 
			// && kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0 
			)
		{
			std::cout << "-"  << " " << rvcKmer
			// << ":" << kmerHash[rvcKmer]->getTotalFreq() 
			<< ":" << kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength);
			// << ":" << kmerHash[rvcKmer]->getTotalSum()/kmerHash[rvcKmer]->getTotalFreq() + m_seedLength
			if(kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength)>0)
				std::cout << ":" << kmerHash[rvcKmer]->getSumOfPos(m_currentLength-m_seedLength)/kmerHash[rvcKmer]->getSumOfFreq(m_currentLength-m_seedLength)+m_seedLength;
		}
		
		if(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end())
		{
			std::cout << "--" << (double)(*iter)->getKmerCount()/m_currentLength;
			std::cout << "\n"; 
		}
		
		// assert(hashiter!=kmerHash.end() || hashiter2!=kmerHash.end());
	}
}

// Print the string represented by every node
void SAIPBSelfCorrectTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node
void SAIPBSelfCorrectTree::attempToExtendUsingHash(STNodePtrList &newLeaves, size_t hashKmerSize, size_t minExtFreq)
{
	double maxLeafFreq = -0.1, removedMaxLeafFreq = -0.1;
	
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		(*iter)->setUpdate(false);
		
		// find the max leaf freq
		double currLeafFreq = (double)(*iter)->getKmerCount()/m_currentLength;
		if(currLeafFreq > maxLeafFreq)
			maxLeafFreq = currLeafFreq;
		
	for(int i=0;i<1;i++)
	{
		std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
		
		// Complex repeats lead to error extensions, requiring at least twice globally.
        extensions = getFMIndexRightExtensions(*iter, minExtFreq);
		// if(extensions.size() <1)
		// {
	    cout<<"extensize:"<<extensions.size()<<"\n";
			// NewextendByNewSeeds(src,(*iter)->getFullString(), hashKmerSize,  maxLength);
			// refineSAInterval(hashKmerSize-1);
			// extensions = getFMIndexRightExtensions(*iter, minExtFreq);
		// }
		
        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
		bool isNoExtension = true;
        if(extensions.size() == 1)
        { 
            // Single extension, do not branch
			// bug:: the leaf may be falsely extended when local kmer freq is insufficient
            // (*iter)->extend(extensions.front().first);
            // (*iter)->fwdInterval=extensions.front().second.interval[0];
            // (*iter)->rvcInterval=extensions.front().second.interval[1];
			
			// std::string fwdkmer = (*iter)->getSuffix(hashKmerSize);
			std::string fwdkmer = (*iter)->getSuffix(hashKmerSize-1)+extensions.front().first;
			
			double currAvgFreq = (double)(*iter)->getKmerCount()/(m_currentLength+1000000);
			
			size_t kmerfreqs = 0;
			if( isExtensionValid(fwdkmer, currAvgFreq, kmerfreqs) )
			{
				(*iter)->setUpdate(true);
				isNoExtension = false;
				(*iter)->extend(extensions.front().first);
				(*iter)->fwdInterval=extensions.front().second.interval[0];
				(*iter)->rvcInterval=extensions.front().second.interval[1];
				(*iter)->addKmerCount(kmerfreqs);
				newLeaves.push_back(*iter);
				
			}
			else if(currLeafFreq > removedMaxLeafFreq)
			
				removedMaxLeafFreq = currLeafFreq;
			
			

        }
        else if(extensions.size() > 1)
        {
            // Create extensions for each new base
			
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                // SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                // pChildNode->fwdInterval=extensions[i].second.interval[0];
                // pChildNode->rvcInterval=extensions[i].second.interval[1];
			
				std::string fwdkmer = (*iter)->getSuffix(hashKmerSize-1)+extensions[i].first;
				
				// std::string fwdkmer = pChildNode->getSuffix(hashKmerSize);
				/**** Math bug: the division requires denominator larger than numerator for accurate average frequency
					numerator denominator division
					100	1000	0.1
					99	999	0.099099099
					100	10	10
					99	9	11
				***/
				double currAvgFreq = (double)(*iter)->getKmerCount()/ (m_currentLength+1000000);
				
				size_t kmerfreqs=0;
				if( isExtensionValid(fwdkmer, currAvgFreq, kmerfreqs) )
				{
					(*iter)->setUpdate(true);
					(*iter)->setUpdate(true);
					isNoExtension = false;
					// inherit kmer freq from parents
					// pChildNode->addKmerCount(kmerfreqs);
					// newLeaves.push_back(pChildNode);
					
					SAIntervalNode* pChildNode= (*iter)->createChild(extensions[i].first);
					pChildNode->fwdInterval=extensions[i].second.interval[0];
					pChildNode->rvcInterval=extensions[i].second.interval[1];
					pChildNode->addKmerCount(kmerfreqs);
					newLeaves.push_back(pChildNode);
				
				}
            }
			
			// find the max leaf freq
			if(isNoExtension && currLeafFreq > removedMaxLeafFreq)
				removedMaxLeafFreq = currLeafFreq;
			
        }
		else	// no extensions
		{
			if(currLeafFreq > removedMaxLeafFreq)
				removedMaxLeafFreq = currLeafFreq;
		}
		// if(!isNoExtension)break;
		// else if (isNoExtension&&i==0)
		// {
			// NewextendByNewSeeds(src,(*iter)->getFullString(), hashKmerSize,  maxLength);
			// refineSAInterval(hashKmerSize-1);
			
		// }
	}
    } // end of each leaf

	if(maxLeafFreq == removedMaxLeafFreq)
		m_isLargeLeaveRemoved = true;

	// std::cout << maxLeafFreq << "\t" << removedMaxLeafFreq <<"\n";
}

// check if new extension satisfies frequency, position, ...
bool SAIPBSelfCorrectTree::isExtensionValid(std::string fwdkmer, double& currAvgFreq, size_t& kmerFreq)
{
		// PB64756_11478, false low-complexity leaves
		// if(isLowComplexity(fwdkmer, 0.8))
			// return false;
			
		kmerHashiter iter1 = kmerHash.find(fwdkmer);			
		
		// bubble removal by removing kmer path with avg kmer freq less than previous one
		// if( iter1!=kmerHash.end() && currAvgFreq < iter1->second->getMaxAvgFreq() ) 
		
		// if(fwdkmer == "GCATCCGGCAA")
			// std::cout << fwdkmer <<": " << currAvgFreq << "\t" << iter1->second->getMaxAvgFreq() << "\n";
		
		// Bubble removal preoduces false removal, perform only when leaves are getting larger than 8
		if( iter1!=kmerHash.end() && m_leaves.size()>8 && currAvgFreq < iter1->second->getMaxAvgFreq() ) 
			return false;
					
		if( iter1!=kmerHash.end() && currAvgFreq > iter1->second->getMaxAvgFreq() )
			iter1->second->setMaxAvgFreq( currAvgFreq );

		// // do it again for reverse complement kmer
		std::string rvckmer = reverseComplement(fwdkmer);
		kmerHashiter iter2 = kmerHash.find(rvckmer);
		
		// Restricted to local kmer frequency
		kmerFreq = iter1==kmerHash.end()? 0 : iter1->second->getSumOfFreq(m_currentLength - m_seedLength);

		// if(fwdkmer == "GCATCCGGCAA")
		// {
			// std::cout << fwdkmer <<": " << kmerFreq << "\t" <<m_currentLength- m_seedLength <<"\n";
			// getchar();
		// }
		
		kmerFreq += iter2==kmerHash.end()? 0 : iter2->second->getSumOfFreq(m_currentLength - m_seedLength);
		
		cout<<kmerFreq<<"kmerFreq\n";
		if(kmerFreq >= m_min_SA_threshold)
			return true;
		else
		{
			cout<<"extenfail!\n";
			return false;
			}
}

// Refine SA intervals of each leave with a new kmer
void SAIPBSelfCorrectTree::refineSAInterval(size_t newKmer)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmer);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));
    }
}

// Convert all leaves back to the status before extension
// trim extended label for single extensions (no children)
// delete children from parents for multiple extensions
void SAIPBSelfCorrectTree::rollBackExtension(size_t newKmerSize, bool btrimAll)
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
		// single extension, trim one base
		if((*iter)->getNumOfChild() == 0 && ((*iter)->isUpdated() || btrimAll) )
		{
			(*iter)->trimLabel(1);
			(*iter)->reduceKmerCount((*iter)->getLastKmerCount());
		}
		else if((*iter)->getNumOfChild() > 0 && ((*iter)->isUpdated() || btrimAll) )
			(*iter)->killAllChild();
			
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmerSize);
        (*iter)->fwdInterval=BWTAlgorithms::findInterval(m_pRBWT, reverse(pkmer));
        (*iter)->rvcInterval=BWTAlgorithms::findInterval(m_pBWT, reverseComplement(pkmer));
    }
}

//  IntervalSizeCutoff min freq at fwd and rvc bwt
std::vector<std::pair<std::string, BWTIntervalPair> > SAIPBSelfCorrectTree::getFMIndexRightExtensions(SAIntervalNode* pNode, const size_t IntervalSizeCutoff)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update forward Interval using extension b
        BWTInterval fwdProbe=pNode->fwdInterval;
        if(fwdProbe.isValid())
            BWTAlgorithms::updateInterval(fwdProbe,b,m_pRBWT);

        //update reverse complement Interval using extension rcb
        BWTInterval rvcProbe=pNode->rvcInterval;
		char rcb=BWT_ALPHABET::getChar(5-i);
        if(rvcProbe.isValid())
            BWTAlgorithms::updateInterval(rvcProbe,rcb,m_pBWT);

        size_t bcount = 0;
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();
		if(bcount >= 1)cout<< b <<"extenword\n";
        if(bcount >= IntervalSizeCutoff)
        {	
			
			// std::cout << m_currentKmerSize << ":" << bcount <<"\n";
            // extend to b
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip;
            bip.interval[0]=fwdProbe;
            bip.interval[1]=rvcProbe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT

    return out;
}


// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIPBSelfCorrectTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;
	STNodePtrList newLeaves;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        assert(currfwd.isValid() || currrvc.isValid());

		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated=currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.lower
                            && currfwd.upper <= m_fwdTerminatedInterval.upper;
        bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.lower
                            && currrvc.upper <= m_rvcTerminatedInterval.upper;

        if(isFwdTerminated || isRvcTerminated)
        {
            std::string STNodeStr = (*iter)->getFullString();
			// std::cout << "Terminated: " << STNodeStr << "\n";
			
            SAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();

            //compute the merged pos right next to the kmer on 2nd read.
            //STresult.index = m_minOverlap ;
            results.push_back(STresult);
            found =  true;
        }
    }

    return found;
}
