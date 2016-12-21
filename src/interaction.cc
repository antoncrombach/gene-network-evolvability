//
// Implementation of a special chromosome element, the interaction.
//
// by Anton Crombach, A.B.M.Crombach at uu.nl
//

#include "interaction.hh"
#include "pool.hh"

template<> blowhole::ObjectCache< blowhole::Interaction >* 
blowhole::ObjectCache< blowhole::Interaction >::instance_ = 0;

blowhole::Interaction::Interaction() 
    : ChromosomeElement(), weight_( -1 ), tag_( 0 ) {}

blowhole::Interaction::Interaction( int t, int w )
    : ChromosomeElement(), weight_( w ), tag_( t ) {}
    
blowhole::Interaction::Interaction( const Interaction &tp ) 
    : ChromosomeElement( tp ) {
    Interaction::copy( tp );
}

blowhole::ChromosomeElement* 
blowhole::Interaction::clone() const {
    //return new Interaction( *this );
    Interaction *ia = ObjectCache< Interaction >::instance()->borrowObject();
    ia->weight_ = weight_;
    ia->tag_ = tag_;
    return ia;
}

void 
blowhole::Interaction::copy( const ChromosomeElement &ce ) {
    const Interaction *ia = dynamic_cast< const Interaction * >( &ce );
    weight_ = ia->weight_;
    tag_ = ia->tag_;
}

bool
blowhole::Interaction::toPool() {
    return ObjectCache< Interaction >::instance()->returnObject( this );
}

void
blowhole::Interaction::writeXml( std::ostream &os ) const { 
    os << "<ia t=\"" << tag_ << "\" w=\"" << weight_ << "\"/>\n";
}

