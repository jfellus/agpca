 /**
 * \file FindFiles.java
 * \author Philippe Gosselin
 * \version 4.0
 */

package retin.toolbox.core;

import java.util.*;

public class BiMap<K,V> {

    HashMap<K,V> keyVal;
    HashMap<V,K> valKey;
    
    public BiMap(){
    }

    public void put(K key, V val){
        if (keyVal == null){
            keyVal = new HashMap<K, V>();
        }
        if (valKey == null){
            valKey = new HashMap<V, K>();
        }

        keyVal.put(key, val);
        valKey.put(val, key);
        
    }

    public K getKey(V val){
        if (valKey == null) return null;
        return valKey.get(val);
    }

    public V getVal(K key){
        if (keyVal == null) return null;
        return keyVal.get(key);
    }

    public boolean containsKey(K key){
        if (keyVal == null) return false;
        return keyVal.containsKey(key);
    }
    
    public boolean containsVal(V val){
        if (valKey == null) return false;
        return valKey.containsKey(val);
    }
}

