//
// XML SAX reader for populations
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_NET_POP_READER_
#define _BLOWHOLE_NET_POP_READER_

#include "defs.hh"
#include "net_agent.hh"
#include "delta_agent.hh"
#include "gene.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace blowhole {

XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    /// \class NetPopReader
    /// \brief Read a population of agents from an xml file and a set of dot
    /// files.
    class NetPopReader : public DefaultHandler {
        public:
        /// ctor
        NetPopReader( Factory * );
        /// dtor
        ~NetPopReader() {};

        /// Return the population
        boost::tuple< std::vector< Agent* >, std::vector< location > >
            population();
        
        /// Read opening tag
        void startElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname,
            const Attributes &attrs );
        /// Read tag contents (not used a lot)
        void characters( const XMLCh* const chars, 
            const unsigned int length );
        /// Read closing tag
        void endElement( const XMLCh* const uri, 
            const XMLCh* const localname, const XMLCh* const qname );
        /// Stuff to zero
        void resetDocument();
    
        /// A few functions for catching warnings and errors
        void error( const SAXParseException &exc );
        void fatalError( const SAXParseException &exc );
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
        Genome *genome_;
        Chromosome *chromo_;
        std::list< ChromosomeElement* > *chr_;
        
        std::vector< Agent* > population_;
        std::vector< location > locations_;
        std::set< int > gene_set_;
        boost::dynamic_bitset<> reference_;
    };

inline bool NetPopReader::sawErrors() const
{ return sawErrors_; }

inline bool NetPopReader::done() const
{ return done_; }

inline void NetPopReader::error( const SAXParseException &e ) 
{ sawErrors_ = true; }

inline void NetPopReader::fatalError( const SAXParseException &e )
{ sawErrors_ = true; }

inline void NetPopReader::resetErrors()
{ sawErrors_ = false; }

inline
boost::tuple< std::vector< blowhole::Agent* >, std::vector< blowhole::location > >
NetPopReader::population() 
{ return boost::make_tuple( population_, locations_ ); }
}
#endif
