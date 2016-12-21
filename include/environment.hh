//
// Environment base class and childs with simple (periodic) changes.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#ifndef _BLOWHOLE_ENVIRONMENT_H_
#define _BLOWHOLE_ENVIRONMENT_H_

#include "defs.hh"
#include "subject.hh"
#include "population.hh"

namespace blowhole {

    /// \class Environment
    /// \brief Abstract base class for modelling the environment.
    ///
    /// The interplay between the environment and the population may become
    /// quite intricate. The current scope of the environment is that it is
    /// some external mechanism influencing the population. The population has
    /// no influence on the environment though.
    class Environment : public Subject {
        public:
        /// Destructor (empty).
        virtual ~Environment() {};
        
        /// Initialise the environment with a random seed
        virtual void initialise( uint );
        /// The environment changes, the population is alerted
        virtual void fluctuate( long ) = 0;
        /// Get attractor to adapt to
        virtual const boost::dynamic_bitset<> & attractor() const = 0;
        /// Finish of the environment
        virtual void finish();
        
        /// Set a reference to the model
        void model( Model * );
        /// Get a reference to the model
        Model* model() const;
        
        /// Has the environment just changed?
        bool hasChanged() const;
        
        protected:
        /// Hidden constructor
        Environment();
        
        protected:
        /// Reference to the model
        Model *model_;
        /// The environment has its own uniform random nr generator
        base_generator_type generator_env_;
        uniform_gen_type uniform_env_;
        /// Flag if environment just changed
        bool changed_;
    };
    
    inline void Environment::model( Model *m )
    { model_ = m; }
    
    inline Model* Environment::model() const
    { return model_; }
    
    inline bool Environment::hasChanged() const
    { return changed_; }
    
    /// \class ConstantEnvironment
    /// \brief Most simple environment, the not-changing one.
    ///
    /// The constant environment is used as a reference to see the influence 
    /// of the other environments.
    class ConstantEnvironment : public Environment {
        public:
        /// Constructor
        ConstantEnvironment( int );
        /// Destructor
        ~ConstantEnvironment() {}
        
        /// No fluctuation in a constant environment.
        virtual void fluctuate( long ) {}
        /// What is optimal...
        virtual const boost::dynamic_bitset<> & attractor() const;
        /// Give attractor of the genes needed to be well adapted to 
        /// this environment
        void attractor( const boost::dynamic_bitset<> & );
        
        private:
        boost::dynamic_bitset<> attractor_;
    };

    inline const boost::dynamic_bitset<> & 
    ConstantEnvironment::attractor() const
    { return attractor_; }

    inline void ConstantEnvironment::attractor( 
        const boost::dynamic_bitset<> &c )
    { attractor_ = c; }

    /// \class PeriodicEnvironment
    /// \brief Environment changes periodically
    ///
    /// Now we know when the environment changes. (Contrast to 
    /// PoissonEnvironment.)
    class PeriodicEnvironment : public Environment {
        public:
        /// Constructor
        PeriodicEnvironment( int );
        /// Destructor
        ~PeriodicEnvironment() {}
        
        /// Initialise every module to the \c low number of gene attractors
        virtual void initialise( uint );
        /// The environment has a chance of switching attractors. 
        void fluctuate( long );
        /// What is optimal...
        virtual const boost::dynamic_bitset<> & attractor() const;
        /// Set the attractors
        void attractor( int, const boost::dynamic_bitset<> & );
        /// Set the length of a period.
        void period( long );
        /// Periods may be shifted to overlap more or less.
        void offset( long );
        
        private:
        uint attractor_index_;
        long period_;
        long offset_;
        
        std::vector< boost::dynamic_bitset<> > states_;
    };

    inline const boost::dynamic_bitset<> &
    PeriodicEnvironment::attractor() const 
    { return states_[ attractor_index_ ]; }

    inline void PeriodicEnvironment::period( long l ) 
    { period_ = l; }

    inline void PeriodicEnvironment::offset( long l ) 
    { offset_ = l; }

    /// \class PoissonEnvironment
    /// \brief Environment changes according to stochastic Poisson process.
    ///
    /// As a short cut to expression levels we take the network attractors as
    /// the criterium for being adapted to the environment. 
    class PoissonEnvironment : public Environment {
        public:
        /// Constructor
        PoissonEnvironment( int );
        /// Destructor
        ~PoissonEnvironment() {}
        
        /// Initialise every module to the \c low number of gene attractors
        virtual void initialise( uint );
        /// Each module had a chance of switching the desired number of
        /// gene attractors. The change is a toggle between two states (values)
        void fluctuate( long );
        /// Get the current `optimal' number of gene attractors
        virtual const boost::dynamic_bitset<> & attractor() const;
        /// Set the attractors
        void attractor( int, const boost::dynamic_bitset<> & );
        /// Set the probability of toggling between two states
        void lambda( double );
        
        private:
        uint attractor_index_;
        double lambda_;
        
        std::vector< boost::dynamic_bitset<> > states_;
    };
    
    inline const boost::dynamic_bitset<> &
    PoissonEnvironment::attractor() const 
    { return states_[ attractor_index_ ]; }

    inline void PoissonEnvironment::lambda( double l ) 
    { lambda_ = l; }
}
#endif
