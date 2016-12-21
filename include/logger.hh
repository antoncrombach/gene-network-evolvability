//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_LOGGER_H_
#define _BLOWHOLE_LOGGER_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace blowhole {

    /// \class LogCsvMutations
    ///
    /// Log the mutations of a timestep. If no mutations, nothing is logged.
    class LogCsvMutations : public AsyncLogObserver {
        public:
        enum clsf { POS, NEG, NEU, UNK };
        
        public:
        LogCsvMutations( std::string );
        virtual ~LogCsvMutations() {}

        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void write();
        void understand( const GenRegAgent & );
        
        private:
        long time_;
        std::vector< uint > cps_, rms_, thrs_, 
            iacps_, iarms_, ians_, iats_, iaws_;
    };    

/*    /// \class LogCsvPerturbations
    ///
    /// Log mean, min, median and stddev of Perturbations.
    class LogCsvPerturbations : public LogObserver {
        public:
        LogCsvPerturbations( std::string, long );
        virtual ~LogCsvPerturbations() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };*/

    /// \class LogGenomeStats
    ///
    /// Gather a bunch of stats and counts from each genome and write to 
    /// some different files different aspects of these genomes.
    ///
    /// How to do this elegantly? Multiple inheritance, this observer
    /// is being observed by three others. Those three are the ones actually
    /// writing to file. This observer notifies them if they need to write...
    class LogGenomeStats : public LogObserver, public Subject {
        public:
        LogGenomeStats( long );
        virtual ~LogGenomeStats() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        void nrGenes( std::vector< double > & ) const;
        void nrActiveGenes( std::vector< double > & ) const;
        void nrRefGenes( std::vector< double > & ) const;
        void nrBindingSites( std::vector< double > & ) const;
        void inDegree( std::vector< double > & ) const;
        void outDegree( std::vector< double > & ) const;
        
        int nrRepeats() const;
        int nrRetroposons() const;
        
        int sample() const;
        
        private:
        void doCounting( Chromosome * );
        void doRefState( const DeltaAgent * );
        
        private:
        int sample_, nr_repeats_, nr_retroposons_;
        std::map< int, int > nr_genes_, nr_genes_on_, nr_bsites_, 
            in_degree_, out_degree_;
        std::vector< int > nr_ref_genes_;
    };

    class LogCsvGenes : public AsyncLogObserver {
        public:
        LogCsvGenes( std::string );
        virtual ~LogCsvGenes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvBindingSites : public AsyncLogObserver {
        public:
        LogCsvBindingSites( std::string );
        virtual ~LogCsvBindingSites() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };
    
    class LogCsvActiveGenes : public AsyncLogObserver {
        public:
        LogCsvActiveGenes( std::string );
        virtual ~LogCsvActiveGenes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvInDegree : public AsyncLogObserver {
        public:
        LogCsvInDegree( std::string );
        virtual ~LogCsvInDegree() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvOutDegree : public AsyncLogObserver {
        public:
        LogCsvOutDegree( std::string );
        virtual ~LogCsvOutDegree() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvRefGenes : public AsyncLogObserver {
        public:
        LogCsvRefGenes( std::string );
        virtual ~LogCsvRefGenes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    class LogCsvRetroposons : public AsyncLogObserver {
        public:
        LogCsvRetroposons( std::string );
        virtual ~LogCsvRetroposons() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };

    /// \class LogCsvDistances
    ///
    /// Log mean, min, median and stddev of distances.
    class LogCsvDistances : public LogObserver {
        public:
        LogCsvDistances( std::string, long );
        virtual ~LogCsvDistances() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvScores
    ///
    /// Log min, mean, median and stddev of fitness scores.
    class LogCsvScores : public LogObserver {
        public:
        LogCsvScores( std::string, long );
        virtual ~LogCsvScores() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvEnvironment
    ///
    /// Log change of environment.
    class LogCsvEnvironment : public AsyncLogObserver {
        public:
        LogCsvEnvironment( std::string );
        virtual ~LogCsvEnvironment() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize() {}
    };        

    /// \class LogCsvAncestors
    ///
    /// Log ancestors to trace them back
    class LogCsvAncestors : public AsyncLogObserver {
        public:
        // ancestor graph
        typedef boost::adjacency_list< 
            boost::listS, boost::setS, boost::bidirectionalS,
            boost::property< boost::vertex_index_t, uint, AgentTag >,
            boost::property< boost::edge_index_t, uint > > ancestor_graph;
        // vertex descriptor
        typedef ancestor_graph::vertex_descriptor anc_vertex;
        // a mapping
        typedef boost::property_map< ancestor_graph, 
            boost::vertex_bundle_t >::type anc_birth_map;
        // accessing vertices
        typedef boost::graph_traits< ancestor_graph >::vertex_iterator
            ver_iter;
        // accessing edges
        typedef boost::graph_traits< ancestor_graph >::adjacency_iterator
            adj_iter;
        typedef boost::inv_adjacency_iterator_generator< ancestor_graph >::type
            invadj_iter;
            
        public:
        LogCsvAncestors( std::string );
        virtual ~LogCsvAncestors() {}
        
        virtual void initialize( Subject * );
        virtual void update( Subject * );
        virtual void finalize();
        
        public:
        static void maxForestSize( uint );
        
        private:
        void writeHeader();
        void writeForest();
        
        void renewForest();
        void initForest();
        
        void birth( Agent * );
        void death( Agent * );
        
        void mapToTime( std::map< AgentTag, anc_vertex > & );
        
        private:
        int nr_roots_;
        anc_vertex root_;
        ancestor_graph forest_;
        std::map< AgentTag, anc_vertex > leafs_, mothers_;
        
        private:
        static uint max_forest_size_;
    };

    /// \class LogXmlAgentTrace
    ///
    /// Log agents given a trace of agents to look for...
    class LogXmlAgentTrace : public AsyncLogObserver {
        // The trace file has a certain format. Unfortunately xml is a bit of
        // a deception, so the trace file is in csv format:
        //
        // each line is an agent id-tag: birth <tab> x <tab> y <nl>
        public:
        LogXmlAgentTrace( std::string, std::string );
        virtual ~LogXmlAgentTrace() {}
        
        virtual void initialize( Subject * );
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        void writeFooter();
        void nextTarget();
        
        private:
        bool first_;
        std::string dname_;
        AgentTag target_;
        boost::filesystem::ifstream *trace_;
    };
    
    /// \class LogXmlGenomes
    ///
    /// Log entire populations.. at least the genomes.
    class LogXmlGenomes : public LogObserver {
        public:
        LogXmlGenomes( std::string, long );
        virtual ~LogXmlGenomes() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    /// \class LogXmlEnvGenomes
    ///
    /// Log entire populations just before the environment is changed.
    class LogXmlEnvGenomes : public AsyncLogObserver {
        public:
        LogXmlEnvGenomes( std::string );
        virtual ~LogXmlEnvGenomes() {}
        
        virtual void initialize( Subject * ) {}
        virtual void update( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    /// \class LogPopulationSize
    ///
    /// Log the number of individuals.
    class LogPopulationSize : public LogObserver {
        public:
        LogPopulationSize( std::string, long );
        virtual ~LogPopulationSize() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };
    
    /// \class LogCsvPopulationDistances
    ///
    /// Log min, mean, median and stddev of distances, 
    /// while being aware of two types of agents.
    class LogCsvPopulationDistances : public LogObserver {
        public:
        LogCsvPopulationDistances( std::string, long );
        virtual ~LogCsvPopulationDistances() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    /// \class LogCsvGrid
    ///
    /// Log the fitness in a matrix format. The resulting file can be used
    /// to generate movies of the grid.
    class LogCsvGrid : public LogObserver {
        public:
        LogCsvGrid( std::string, long );
        virtual ~LogCsvGrid() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

}
#endif

