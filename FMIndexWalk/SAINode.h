//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

// 
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAINODE_H
#define SAINODE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"

class SAINode;
typedef std::list<SAINode*> SAINodePtrList;

// Base class represent search nodes for all sorts of applications using FM-index walk
class SAINode
{
    public:

        //
        // Functions
        //
        SAINode(const std::string* pQuery, SAINode* parent);
        ~SAINode();

        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAINode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);

        // Return a suffix of length l of the string represented by this node
        std::string getSuffix(size_t l) const;

        // Return the complete sequence of the string represented by the branch
        std::string getFullString() const;

        // Initialize or update the alignment data
        //void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeInitial(const std::string& initialLabel);

        // Recursive function to print all the strings represented
        // by this node and all its children.
        void printAllStrings(const std::string& parent) const;

        size_t getKmerCount(){return m_totalKmerCount;};
		
		size_t getLastKmerCount(){return m_lastKmerCount;};
		
        void addKmerCount(size_t currKmerCount){ 
			m_totalKmerCount += currKmerCount; 
			m_lastKmerCount = currKmerCount;
		};

        void reduceKmerCount(size_t currKmerCount){ 
			m_totalKmerCount -= currKmerCount; 
			// m_lastKmerCount = currKmerCount;
		};
		
		// get number of children
		size_t getNumOfChild(){return m_children.size();}
		
		// kill all children
		void killAllChild()
		{
			m_children.clear();
		}
		
		// trim label
		void trimLabel(size_t trimLen)
		{
			m_label = m_label.substr(0, m_label.length()-trimLen);
		}
		
		inline void setUpdate(bool flag)
		{
			m_updated = flag;
		}
		
		inline bool isUpdated(){return m_updated;};

    private:

        //
        // Data
        //

        // The extension string from the parent
        std::string m_label;
        size_t m_totalKmerCount;
		size_t m_lastKmerCount;
		bool m_updated;

        // size_t m_kmersize;


	protected:
        // The query string being threaded through the graph
        const std::string* m_pQuery;

        // The parent node, can be NULL
        SAINode* m_pParent;
        SAINodePtrList m_children;

};


//
// SAIntervalNode for implementation of FM-index walk between paired end reads or validation
//
class SAIntervalNode : public SAINode
{
    public:

        //
        // Functions
        //
        SAIntervalNode(const std::string* pQuery, SAIntervalNode* parent):SAINode(pQuery,parent)
		{
		}
        ~SAIntervalNode(){};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIntervalNode* createChild(const std::string& label);

        BWTInterval fwdInterval;
        BWTInterval rvcInterval;

};

// leaves of SAIntervalNode
typedef std::list<SAIntervalNode*> STNodePtrList;

// Object to hold the result of the extension process
struct SAIntervalNodeResult
{
    std::string thread;
	size_t SAICoverage;
	int SAIntervalSize;
};
typedef std::vector<SAIntervalNodeResult> SAIntervalNodeResultVector;



//
// SAIOverlapNode for implementation of overlap computation using FM-index walk
//
class SAIOverlapNode : public SAINode
{
    public:

        //
        // Functions
        //
        SAIOverlapNode(const std::string* pQuery, SAIOverlapNode* parent):SAINode(pQuery,parent)
		{
			lastSeedIdx=totalSeeds=lastOverlapLen=currOverlapLen=queryOverlapLen=numOfErrors=0;
			lastSeedIdxOffset=0;
		}
        ~SAIOverlapNode(){};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIOverlapNode* createChild(const std::string& label);

        BWTIntervalPair currIntervalPair;
		
		// last matched seed index
		size_t lastSeedIdx;
		// last overlap length when matching last seed
		size_t lastOverlapLen;
		size_t totalSeeds;
		// current overlap length on the subject increases wrt each FM-index extension
		size_t currOverlapLen;
		// number of SNPs or indels
		size_t numOfErrors;
		// index offset to the center
		int lastSeedIdxOffset;

		// index of the init seed
		int initSeedIdx;
		
		// current overlap length on the query
		size_t queryOverlapLen;

};

// leaves of SAIOverlapNode
typedef std::list<SAIOverlapNode*> SONodePtrList;

#endif
