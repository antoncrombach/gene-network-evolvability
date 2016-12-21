//
// XML SAX reader for agents
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_NET_AGENT_READER_
#define _BLOWHOLE_NET_AGENT_READER_

#include "defs.hh"
#include "net_agent.hh"
#include "delta_agent.hh"
#include "gene.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace blowhole {

// Xerces shit.
XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    /// \class NetAgentReader
    /// \brief Read in an agent from an xml and dot file.
    ///
    /// NetAgentReader reads a genome from an xml file and a network from a dot
    /// file.
    class NetAgentReader : public DefaultHandler {
        public:
        /// Constructor.
        NetAgentReader( Factory *, int );
        /// Destructor.
        ~NetAgentReader() {};

        // Return the agent
        GenRegAgent* agent();
        
        /// Read opening xml tag.
        void startElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname,
            const Attributes &attrs );
        /// Read tag contents. 
        void characters( const XMLCh* const chars, const uint length );
        /// Read closing xml tag.
        void endElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname );
        /// Set stuff to zero 
        void resetDocument();
    
        /// A few functions for catching warnings and errors
        void error( const SAXParseException & );
        void fatalError( const SAXParseException & );
        void resetErrors();
        bool sawErrors() const;
        bool done() const;
        
        private:
        Factory *factory_;
        bool sawErrors_, done_, class_, ref_state_;
        
        std::string agent_class_;
        int type_;
        
        GenRegAgent *agent_;
        AgentTag me_, mother_;
        HopfieldGraph *network_;
        Genome *genome_;
        Chromosome *chromo_;
        std::list< ChromosomeElement* > *chr_;
        std::set< int > gene_set_;
        boost::dynamic_bitset<> reference_;
    };

inline bool NetAgentReader::sawErrors() const
{ return sawErrors_; }

inline bool NetAgentReader::done() const
{ return done_; }

inline void NetAgentReader::error( const SAXParseException &e ) 
{ sawErrors_ = true; }

inline void NetAgentReader::fatalError( const SAXParseException &e )
{ sawErrors_ = true; }

inline void NetAgentReader::resetErrors()
{ sawErrors_ = false; }

inline GenRegAgent* NetAgentReader::agent() 
{ return agent_; }
}
#endif
